#!/usr/bin/env python3
"""Chain co-occurrence graph for the four-donor Immudex dCODE matrices.

A node is one chain (locus + V + J + CDR3). An edge is two chains in the same
cell barcode, carrying both UMI counts. Leakage is a Poisson ambient component

    lambda(b, c) = rho_b * U_c / sum(U)

with rho_b fit once from the chains in barcode b that cannot be the primary
receptor (a third TRA, a second TRB). A chain is endogenous when that mixture
assigns it to the high component. Kept clones are alpha-beta or alpha-alpha-beta.

The same graph names the parent of a leaked chain. A fraction of that parent's
feature profile is subtracted from the barcode. A clone is then specific for a
pMHC when its mean per-cell residual exceeds every per-cell negative-control
count from that donor and the peptide's clone-mean histogram has a local mode
at or above that count. A sticky mode whose tail only leaks past the controls
is not a call. Clone size cannot sum a mild residual over the bar.
There is no UMI cutoff.

The on-disk "binarized" matrix is numeric counts. Columns ending in ``_binder``
are a previous 0/1 call and are not read. There is no transcriptome in the
dCODE directory; surface-protein counts stand in for the expression background.

# 2026-10-16. Issue https://github.com/antigenomics/vdjdb-db/issues/1357
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import polars as pl

CELL_FLAGS = ("True", "true", "TRUE")
PRODUCTIVE = ("True", "true", "TRUE")
LOCI = ("TRA", "TRB")
LAM_FLOOR = 1e-8
META_COLUMNS = ("barcode", "donor", "cell_clono_cdr3_aa", "cell_clono_cdr3_nt")
ALLELE_RE = re.compile(r"^([ABC])(\d{2})(\d{2})$")


def chain_table(contigs: pl.DataFrame, donor: int) -> pl.DataFrame:
    """One row per chain per barcode. Productive TRA/TRB only.

    ``is_cell`` false barcodes stay in the table so their UMIs enter the
    ambient frequency. They are not clone barcodes. ``high_confidence`` is not
    a filter: low-confidence contigs are the leakage observations.
    """
    if contigs.schema["productive"] == pl.Boolean:
        productive_ok = pl.col("productive").fill_null(False)
    else:
        productive_ok = pl.col("productive").cast(pl.String).is_in(list(PRODUCTIVE))
    df = contigs.filter(
        pl.col("chain").is_in(list(LOCI))
        & productive_ok
        & pl.col("cdr3").is_not_null()
        & (pl.col("cdr3") != "")
        & (pl.col("cdr3") != "None")
    )
    if df.schema["is_cell"] == pl.Boolean:
        cell_flag = pl.col("is_cell")
    else:
        cell_flag = pl.col("is_cell").cast(pl.String).is_in(list(CELL_FLAGS))
    df = df.with_columns(
        pl.lit(donor).alias("donor"),
        pl.col("umis").cast(pl.Float64).alias("umi"),
        cell_flag.alias("is_cell"),
        (
            pl.col("chain")
            + pl.lit("|")
            + pl.col("v_gene").fill_null("")
            + pl.lit("|")
            + pl.col("j_gene").fill_null("")
            + pl.lit("|")
            + pl.col("cdr3")
        ).alias("chain_id"),
        pl.col("chain").alias("locus"),
        pl.col("v_gene").alias("v_call"),
        pl.col("j_gene").alias("j_call"),
        pl.col("cdr3").alias("cdr3"),
    )
    return (
        df.group_by("donor", "barcode", "chain_id", "locus", "v_call", "j_call", "cdr3")
        .agg(pl.col("umi").sum(), pl.col("is_cell").max())
    )


def edge_table(chains: pl.DataFrame) -> pl.DataFrame:
    """Undirected co-occurrence edges inside cell barcodes."""
    cell = chains.filter(pl.col("is_cell")).select(
        "donor",
        "barcode",
        pl.col("chain_id").alias("chain_a"),
        pl.col("locus").alias("locus_a"),
        pl.col("umi").alias("umi_a"),
    )
    other = chains.filter(pl.col("is_cell")).select(
        "donor",
        "barcode",
        pl.col("chain_id").alias("chain_b"),
        pl.col("locus").alias("locus_b"),
        pl.col("umi").alias("umi_b"),
    )
    return (
        cell.join(other, on=["donor", "barcode"])
        .filter(pl.col("chain_a") < pl.col("chain_b"))
    )


def _fit_excess(y: np.ndarray, lam: np.ndarray, n_iter: int = 25) -> tuple[float, float]:
    """Moment estimate of the endogenous excess mu and its mixing weight pi.

    Ambient mean lam is known per barcode. The endogenous component is Poisson
    with mean lam + mu. Responsibilities use that likelihood; mu is the weighted
    mean of (y - lam) on the endogenous component, floored at 0. One fit, no
    cutoff grid.
    """
    y = np.asarray(y, dtype=np.float64)
    lam = np.maximum(np.asarray(lam, dtype=np.float64), LAM_FLOOR)
    gap = y - lam
    mu = max(float(np.quantile(gap, 0.9)), 1.0)
    pi = 0.2
    for _ in range(n_iter):
        rate_e = lam + mu
        log_e = math.log(pi) + y * np.log(rate_e) - rate_e
        log_a = math.log(1.0 - pi) + y * np.log(lam) - lam
        # log y! cancels. Shift for stability.
        shift = np.maximum(log_e, log_a)
        re = np.exp(log_e - shift)
        ra = np.exp(log_a - shift)
        r = re / (re + ra)
        pi = float(np.clip(r.mean(), 1e-3, 1.0 - 1e-3))
        weight = float(r.sum())
        if weight < 1e-6:
            return 0.0, pi
        mu_next = max(float((r * gap).sum() / weight), 0.0)
        if abs(mu_next - mu) < 1e-3:
            mu = mu_next
            break
        mu = mu_next
    return mu, pi


def mark_endogenous(chains: pl.DataFrame) -> pl.DataFrame:
    """Add ``endogenous`` to each cell-barcode chain row.

    Non-cell rows are returned with ``endogenous`` false. They have already
    contributed to the chain frequency f_c.
    """
    if chains.is_empty():
        return chains.with_columns(pl.lit(False).alias("endogenous"))

    totals = chains.group_by("chain_id").agg(pl.col("umi").sum().alias("u_c"))
    u_sum = float(totals.select(pl.col("u_c").sum()).item())
    totals = totals.with_columns((pl.col("u_c") / u_sum).alias("f_c"))

    cells = chains.filter(pl.col("is_cell"))
    if cells.is_empty():
        return chains.with_columns(pl.lit(False).alias("endogenous"))

    # Ties in UMI break by chain id. Ordinal rank on an unsorted frame follows
    # row order, and that order is not stable across runs.
    ranked = cells.sort(
        ["barcode", "locus", "umi", "chain_id"],
        descending=[False, False, True, False],
    ).with_row_index("_ord")
    first = ranked.group_by("barcode", "locus").agg(pl.col("_ord").min().alias("_ord0"))
    ranked = ranked.join(first, on=["barcode", "locus"]).with_columns(
        (pl.col("_ord") - pl.col("_ord0") + 1).alias("locus_rank")
    ).drop("_ord", "_ord0").with_columns(
        (
            ((pl.col("locus") == "TRA") & (pl.col("locus_rank") <= 2))
            | ((pl.col("locus") == "TRB") & (pl.col("locus_rank") <= 1))
        ).alias("primary_slot"),
        (
            ((pl.col("locus") == "TRA") & (pl.col("locus_rank") > 2))
            | ((pl.col("locus") == "TRB") & (pl.col("locus_rank") > 1))
        ).alias("excess_slot"),
    )
    ranked = ranked.join(totals.select("chain_id", "f_c", "u_c"), on="chain_id", how="left")
    # A chain seen in one barcode is not the soup. Using it as rho would force
    # lam equal to its own count and hide a second receptor.
    n_per = ranked.group_by("chain_id").agg(pl.len().alias("n_cell_all"))
    ranked = ranked.join(n_per, on="chain_id", how="left")

    excess = (
        ranked.filter(pl.col("excess_slot") & (pl.col("n_cell_all") >= 2))
        .group_by("barcode")
        .agg(pl.col("umi").sum().alias("excess_umi"), pl.col("f_c").sum().alias("excess_f"))
        .with_columns(
            pl.when(pl.col("excess_f") > 0)
            .then(pl.col("excess_umi") / pl.col("excess_f"))
            .otherwise(None)
            .alias("rho")
        )
    )
    rho_fill = excess.select(pl.col("rho").drop_nulls().median()).item()
    if rho_fill is None:
        rho_fill = 0.0
    ranked = ranked.join(excess.select("barcode", "rho"), on="barcode", how="left").with_columns(
        pl.col("rho").fill_null(float(rho_fill)),
    )
    ranked = ranked.with_columns((pl.col("rho") * pl.col("f_c")).clip(lower_bound=LAM_FLOOR).alias("lam"))

    # Abundant chains get one mixture. Rare chains use their own primary excess:
    # a chain seen in a handful of barcodes cannot be the ambient soup.
    sizes = ranked.group_by("chain_id").agg(pl.len().alias("n_cell"), pl.col("u_c").first())
    abundant_ids = sizes.filter((pl.col("n_cell") >= 8) & (pl.col("u_c") >= 20))["chain_id"].to_list()
    mu_pi: dict[str, tuple[float, float]] = {}
    if abundant_ids:
        sub = ranked.filter(pl.col("chain_id").is_in(abundant_ids)).sort("chain_id")
        ids = sub["chain_id"].to_numpy()
        y = sub["umi"].to_numpy()
        lam = sub["lam"].to_numpy()
        if len(ids):
            cuts = np.flatnonzero(ids[1:] != ids[:-1]) + 1
            starts = np.concatenate(([0], cuts))
            ends = np.concatenate((cuts, [len(ids)]))
            for start, end in zip(starts, ends, strict=True):
                mu, pi = _fit_excess(y[start:end], lam[start:end])
                mu_pi[str(ids[start])] = (mu, pi)

    primary = ranked.filter(pl.col("primary_slot"))
    own = primary.group_by("chain_id").agg(
        pl.len().alias("n_primary"),
        (pl.col("umi") - pl.col("lam")).max().alias("mu_max"),
    )
    n_cell = ranked.group_by("chain_id").agg(pl.len().alias("n_cell"))
    depth = ranked.group_by("barcode").agg(pl.col("umi").max().alias("depth"))
    med_depth = float(depth.select(pl.col("depth").median()).item() or 1.0)
    pos = primary.filter((pl.col("umi") - pl.col("lam")) > 0)
    mu_global = 1.0 if pos.is_empty() else float(
        pos.select((pl.col("umi") - pl.col("lam")).median()).item()
    )
    depth = depth.with_columns((pl.col("depth") / med_depth).alias("depth_scale"))

    fitted = pl.DataFrame(
        {
            "chain_id": list(mu_pi),
            "mu_fit": [v[0] for v in mu_pi.values()],
            "pi_fit": [v[1] for v in mu_pi.values()],
        }
    ) if mu_pi else pl.DataFrame(
        schema={"chain_id": pl.String, "mu_fit": pl.Float64, "pi_fit": pl.Float64}
    )

    ranked = (
        ranked.join(own, on="chain_id", how="left")
        .join(n_cell, on="chain_id", how="left")
        .join(depth.select("barcode", "depth_scale"), on="barcode", how="left")
        .join(fitted, on="chain_id", how="left")
        .with_columns(
            pl.col("n_primary").fill_null(0),
            pl.col("mu_max").fill_null(0.0),
        )
    )
    # π pseudocount keeps a chain that is primary in every barcode able to lose
    # a low-UMI member, and keeps an unseen chain from a log(0).
    ranked = ranked.with_columns(
        pl.when(pl.col("mu_fit").is_not_null())
        .then(pl.col("mu_fit"))
        .when(pl.col("n_primary") >= 1)
        .then(pl.col("mu_max").clip(lower_bound=0.0))
        .otherwise(mu_global * pl.col("depth_scale"))
        .alias("mu"),
        pl.when(pl.col("pi_fit").is_not_null())
        .then(pl.col("pi_fit"))
        .otherwise((pl.col("n_primary") + 0.5) / (pl.col("n_cell") + 1.0))
        .alias("pi"),
    )
    ranked = ranked.with_columns(
        (
            (pl.col("pi") / (1.0 - pl.col("pi"))).log()
            + pl.col("umi") * ((pl.col("lam") + pl.col("mu").clip(lower_bound=1e-6)) / pl.col("lam")).log()
            - pl.col("mu")
        ).alias("log_odds")
    )
    # A count that does not clear the ambient mean is the ambient component.
    ranked = ranked.with_columns(
        ((pl.col("log_odds") > 0) & (pl.col("umi") > pl.col("lam"))).alias("endogenous")
    )
    keep = [
        "donor", "barcode", "chain_id", "locus", "v_call", "j_call", "cdr3",
        "umi", "is_cell", "rho", "lam", "mu", "pi", "log_odds", "endogenous",
        "locus_rank",
    ]
    return ranked.select(keep)


def pairing_label(n_tra: int, n_trb: int) -> str:
    """scirpy ``chain_pairing`` names, from the productive chain counts.

    Recorded so the kept dual-alpha set can be compared with ``extra VJ``. Not a filter.
    """
    if n_tra > 2 or n_trb > 2:
        return "multichain"
    if n_tra == 2 and n_trb == 2:
        return "two full chains"
    if n_tra == 2 and n_trb == 1:
        return "extra VJ"
    if n_tra == 1 and n_trb == 2:
        return "extra VDJ"
    if n_tra == 1 and n_trb == 1:
        return "single pair"
    if n_tra == 1 and n_trb == 0:
        return "orphan VJ"
    if n_tra == 0 and n_trb == 1:
        return "orphan VDJ"
    if n_tra == 2 and n_trb == 0:
        return "orphan VJ"
    if n_tra == 0 and n_trb == 2:
        return "orphan VDJ"
    return "no IR"


def clone_table(endogenous: pl.DataFrame) -> pl.DataFrame:
    """Barcodes whose endogenous chains are exactly alpha-beta or alpha-alpha-beta."""
    keep = endogenous.filter(pl.col("endogenous"))
    if keep.is_empty():
        return pl.DataFrame(
            schema={
                "donor": pl.Int64,
                "barcode": pl.String,
                "shape": pl.String,
                "alpha1": pl.String,
                "alpha2": pl.String,
                "beta": pl.String,
                "clone_id": pl.String,
            }
        )
    tra = (
        keep.filter(pl.col("locus") == "TRA")
        .sort("chain_id")
        .group_by("donor", "barcode")
        .agg(pl.col("chain_id").alias("alphas"), pl.len().alias("n_tra"))
    )
    trb = (
        keep.filter(pl.col("locus") == "TRB")
        .group_by("donor", "barcode")
        .agg(pl.col("chain_id").alias("betas"), pl.len().alias("n_trb"))
    )
    both = tra.join(trb, on=["donor", "barcode"], how="inner")
    both = both.filter(
        ((pl.col("n_tra") == 1) & (pl.col("n_trb") == 1))
        | ((pl.col("n_tra") == 2) & (pl.col("n_trb") == 1))
    )
    both = both.with_columns(
        pl.when(pl.col("n_tra") == 1).then(pl.lit("ab")).otherwise(pl.lit("aab")).alias("shape"),
        pl.col("alphas").list.get(0).alias("alpha1"),
        pl.when(pl.col("n_tra") == 2).then(pl.col("alphas").list.get(1)).otherwise(None).alias("alpha2"),
        pl.col("betas").list.get(0).alias("beta"),
    )
    both = both.with_columns(
        pl.when(pl.col("shape") == "ab")
        .then(pl.col("alpha1") + pl.lit("||") + pl.col("beta"))
        .otherwise(pl.col("alpha1") + pl.lit("|") + pl.col("alpha2") + pl.lit("||") + pl.col("beta"))
        .alias("clone_id")
    )
    return both.select("donor", "barcode", "shape", "alpha1", "alpha2", "beta", "clone_id")


def _feature_matrix(features: pl.DataFrame, columns: list[str]) -> tuple[np.ndarray, dict[str, int]]:
    barcodes = features["barcode"].to_list()
    index = {barcode: i for i, barcode in enumerate(barcodes)}
    matrix = features.select(columns).to_numpy().astype(np.float64)
    return matrix, index


def _poisson_ll(y: np.ndarray, rate: np.ndarray) -> float:
    rate = np.maximum(rate, 1e-3)
    return float(np.sum(y * np.log(rate) - rate))


def drop_doublets(
    clones: pl.DataFrame,
    features: pl.DataFrame,
    count_cols: list[str],
) -> pl.DataFrame:
    """Drop an alpha-alpha-beta barcode when one alpha has another beta and the counts fit the sum.

    Parents are the alpha-beta clones of those chains, not random cells. A second
    alpha that never forms another clone stays: that is a dual-alpha receptor.
    """
    if clones.is_empty() or "aab" not in clones["shape"].to_list():
        return clones
    # Barcode order makes the parent-mean sums and the partner lists reproducible.
    ab = clones.filter(pl.col("shape") == "ab").sort("barcode")
    partners: dict[str, list[tuple[str, str]]] = {}
    for alpha, beta, clone_id in zip(
        ab["alpha1"].to_list(), ab["beta"].to_list(), ab["clone_id"].to_list(), strict=True
    ):
        partners.setdefault(alpha, []).append((beta, clone_id))
    # One profile per alpha-beta clone.
    matrix, index = _feature_matrix(features, count_cols)
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, int] = {}
    for barcode, clone_id in zip(ab["barcode"].to_list(), ab["clone_id"].to_list(), strict=True):
        row = index.get(barcode)
        if row is None:
            continue
        if clone_id not in sums:
            sums[clone_id] = np.zeros(len(count_cols), dtype=np.float64)
            counts[clone_id] = 0
        sums[clone_id] += matrix[row]
        counts[clone_id] += 1
    means = {clone_id: sums[clone_id] / counts[clone_id] for clone_id in sums}

    drop: set[str] = set()
    aab = clones.filter(pl.col("shape") == "aab").sort("barcode")
    for barcode, a1, a2, beta in zip(
        aab["barcode"].to_list(),
        aab["alpha1"].to_list(),
        aab["alpha2"].to_list(),
        aab["beta"].to_list(),
        strict=True,
    ):
        foreign = []
        for alpha in (a1, a2):
            for other_beta, clone_id in partners.get(alpha, []):
                if other_beta != beta and clone_id in means:
                    foreign.append(clone_id)
        if not foreign:
            continue
        # Both alphas already head a different beta clone: two receptors in one droplet.
        a1_foreign = any(b != beta for b, _cid in partners.get(a1, []))
        a2_foreign = any(b != beta for b, _cid in partners.get(a2, []))
        row_i = index.get(barcode)
        if row_i is None or (a1_foreign and a2_foreign):
            drop.add(barcode)
            continue
        parent_ids = sorted(set(foreign))
        y = matrix[row_i]
        # The alpha-beta clone of this barcode's own beta, if one alpha also forms it.
        own_ids = sorted({
            clone_id
            for alpha in (a1, a2)
            for other_beta, clone_id in partners.get(alpha, [])
            if other_beta == beta and clone_id in means
        })
        if not own_ids:
            # No separate alpha-beta profile for this pair. Drop when the sum of two
            # foreign parents fits better than either parent alone.
            for i, left in enumerate(parent_ids):
                for right in parent_ids[i + 1:]:
                    both = means[left] + means[right]
                    if (
                        _poisson_ll(y, both) > _poisson_ll(y, means[left])
                        and _poisson_ll(y, both) > _poisson_ll(y, means[right])
                    ):
                        drop.add(barcode)
            continue
        for own_id in own_ids:
            own = means[own_id]
            for other_id in parent_ids:
                if _poisson_ll(y, own + means[other_id]) > _poisson_ll(y, own):
                    drop.add(barcode)
                    break
    if not drop:
        return clones
    return clones.filter(~pl.col("barcode").is_in(list(drop)))


def subtract_leakage(
    endogenous: pl.DataFrame,
    clones: pl.DataFrame,
    features: pl.DataFrame,
    count_cols: list[str],
) -> pl.DataFrame:
    """Subtract each leaked chain's parent profile from kept barcodes.

    Fraction is this barcode's UMI of the leaked chain divided by that chain's
    mean endogenous UMI. The parent is the kept clone in which the chain's
    endogenous UMI is highest.
    """
    corrected = features.select(["barcode", *count_cols]).with_columns(
        [pl.col(c).cast(pl.Float64) for c in count_cols]
    )
    if clones.is_empty():
        return corrected
    kept = set(clones["barcode"].to_list())
    members = clones.select(
        "barcode", "clone_id", "alpha1", "alpha2", "beta"
    )
    long = (
        members.unpivot(
            index=["barcode", "clone_id"],
            on=["alpha1", "alpha2", "beta"],
            value_name="chain_id",
        )
        .drop("variable")
        .drop_nulls("chain_id")
    )
    endo_umi = endogenous.filter(pl.col("endogenous")).select("barcode", "chain_id", "umi")
    long = long.join(endo_umi, on=["barcode", "chain_id"], how="left")
    chain_parent = (
        long.group_by("chain_id", "clone_id")
        .agg(pl.col("umi").mean().alias("mean_umi"), pl.len().alias("n"))
        .sort(["mean_umi", "clone_id"], descending=[True, False])
        .unique(subset=["chain_id"], keep="first")
    )
    parent_of = dict(zip(chain_parent["chain_id"].to_list(), chain_parent["clone_id"].to_list(), strict=True))
    mean_umi = dict(zip(chain_parent["chain_id"].to_list(), chain_parent["mean_umi"].to_list(), strict=True))

    matrix, index = _feature_matrix(corrected, count_cols)
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, int] = {}
    ordered = clones.sort("barcode")
    for barcode, clone_id in zip(ordered["barcode"].to_list(), ordered["clone_id"].to_list(), strict=True):
        row = index.get(barcode)
        if row is None:
            continue
        if clone_id not in sums:
            sums[clone_id] = np.zeros(len(count_cols))
            counts[clone_id] = 0
        sums[clone_id] += matrix[row]
        counts[clone_id] += 1
    means = {cid: sums[cid] / counts[cid] for cid in sums}

    present = (
        endogenous.filter(pl.col("barcode").is_in(list(kept)) & ~pl.col("endogenous"))
        .sort("barcode", "chain_id")
    )
    for barcode, chain_id, umi in zip(
        present["barcode"].to_list(), present["chain_id"].to_list(), present["umi"].to_list(), strict=True
    ):
        parent = parent_of.get(chain_id)
        row = index.get(barcode)
        base = mean_umi.get(chain_id) or 0.0
        if parent is None or row is None or base <= 0 or parent not in means:
            continue
        fraction = float(umi) / float(base)
        matrix[row] = np.maximum(matrix[row] - fraction * means[parent], 0.0)

    out = pl.DataFrame({"barcode": corrected["barcode"].to_list()})
    for j, col in enumerate(count_cols):
        out = out.with_columns(pl.Series(col, matrix[:, j]))
    return out


def negative_control_columns(columns: list[str]) -> list[str]:
    return [c for c in columns if "_NC" in c and not c.endswith("_binder")]


def dextramer_columns(columns: list[str]) -> list[str]:
    out = []
    for col in columns:
        if col in META_COLUMNS or col.endswith("_binder") or "_NC" in col:
            continue
        head = col.split("_", 1)[0]
        if ALLELE_RE.match(head):
            out.append(col)
    return out


def surface_columns(columns: list[str]) -> list[str]:
    skip = set(META_COLUMNS)
    dex = set(dextramer_columns(columns))
    nc = set(negative_control_columns(columns))
    return [
        c for c in columns
        if c not in skip and not c.endswith("_binder") and c not in dex and c not in nc
    ]


def column_residuals(
    corrected: pl.DataFrame, dex_cols: list[str], nc_cols: list[str],
) -> tuple[pl.DataFrame, dict[str, float]]:
    """Residual dextramer counts after an additive and proportional control effect.

    alpha_p is the median count where the negative-control total is zero: the
    floor a sticky dextramer has even with no control signal. gamma_p is the
    median of (count - alpha) / control total where the control total is
    positive. A rare peptide moves neither median.
    """
    if not nc_cols:
        residuals = corrected.select(["barcode", *dex_cols])
        return residuals, {col: 0.0 for col in dex_cols}
    work = corrected.with_columns(
        pl.sum_horizontal([pl.col(c) for c in nc_cols]).alias("nc_total")
    )
    silent = work.filter(pl.col("nc_total") == 0)
    loaded = work.filter(pl.col("nc_total") > 0)
    gammas: dict[str, float] = {}
    alphas: dict[str, float] = {}
    for col in dex_cols:
        alpha = 0.0 if silent.is_empty() else float(silent.select(pl.col(col).median()).item() or 0.0)
        if loaded.is_empty():
            gamma = 0.0
        else:
            gamma = loaded.select(
                ((pl.col(col) - alpha) / pl.col("nc_total")).median()
            ).item()
            gamma = max(float(gamma or 0.0), 0.0)
        alphas[col] = alpha
        gammas[col] = gamma
    residual = work.select("barcode")
    for col in dex_cols:
        expected = alphas[col] + gammas[col] * work["nc_total"]
        residual = residual.with_columns((work[col] - expected).clip(lower_bound=0.0).alias(col))
    return residual, gammas


def _empty_calls() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "donor": pl.Int64,
            "clone_id": pl.String,
            "n_cells": pl.Int64,
            "column": pl.String,
            "peptide": pl.String,
            "mhc": pl.String,
            "residual": pl.Float64,
        }
    )


def control_envelope(clones: pl.DataFrame, controls: pl.DataFrame, nc_cols: list[str]) -> float:
    """Maximum per-cell negative-control count on the kept clone barcodes.

    A specificity has to fall outside the control sample, not above a quantile
    of it. A sticky control reagent sets the bar. This is not a scanned UMI cutoff.
    """
    if clones.is_empty() or not nc_cols:
        return 0.0
    joined = clones.select("barcode").unique().join(
        controls.select(["barcode", *nc_cols]), on="barcode", how="inner",
    )
    if joined.is_empty():
        return 0.0
    ceiling = joined.select(pl.max_horizontal(nc_cols).max()).item()
    return float(ceiling or 0.0)


def component_threshold(values: np.ndarray, envelope: float) -> float | None:
    """Call threshold for one peptide, or None when the tail is not its own mode.

    Bins are powers of two. A local mode is a bin above both neighbors. Only a
    mode whose left edge is at or above the control envelope qualifies; the
    threshold is that edge, walked left along the ascent while the bin stays
    above the envelope. A mode centered below the envelope is stickiness, and
    its leak past the controls is not called.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return None
    top = float(max(values.max(), envelope, 1.0))
    edges = [0.0, 1.0]
    while edges[-1] <= top:
        edges.append(edges[-1] * 2.0)
    hist, _ = np.histogram(values, bins=edges)
    modes = []
    for i in range(1, len(hist)):
        right = hist[i + 1] if i + 1 < len(hist) else -1
        if hist[i] > hist[i - 1] and hist[i] >= right:
            modes.append(i)
    qualify = [i for i in modes if edges[i] >= envelope]
    if not qualify:
        return None
    start = qualify[0]
    while start > 0 and hist[start - 1] < hist[start] and edges[start - 1] >= envelope:
        start -= 1
    return max(envelope, edges[start])


def call_specificities(
    clones: pl.DataFrame,
    residuals: pl.DataFrame,
    dex_cols: list[str],
    controls: pl.DataFrame,
    nc_cols: list[str],
) -> pl.DataFrame:
    """One row per clone per pMHC whose mean residual clears the control envelope and its mode."""
    if clones.is_empty():
        return _empty_calls()
    ceiling = control_envelope(clones, controls, nc_cols)
    joined = clones.select("donor", "barcode", "clone_id").join(residuals, on="barcode", how="inner")
    agg = joined.group_by("donor", "clone_id").agg(
        pl.len().alias("n_cells"),
        *[pl.col(c).mean().alias(c) for c in dex_cols],
    )
    rows: list[dict[str, object]] = []
    donor = agg["donor"].to_numpy()
    clone_ids = agg["clone_id"].to_list()
    n_cells = agg["n_cells"].to_numpy()
    for col in dex_cols:
        peptide, mhc = parse_dextramer(col)
        values = agg[col].to_numpy()
        threshold = component_threshold(values, ceiling)
        if threshold is None:
            continue
        for i in np.flatnonzero(values > threshold):
            rows.append({
                "donor": int(donor[i]),
                "clone_id": clone_ids[i],
                "n_cells": int(n_cells[i]),
                "column": col,
                "peptide": peptide,
                "mhc": mhc,
                "residual": float(values[i]),
            })
    if not rows:
        return _empty_calls()
    return pl.DataFrame(rows)


def parse_dextramer(column: str) -> tuple[str, str]:
    """``A0301_KLGGALQAK_IE-1_CMV`` → peptide, ``HLA-A*03:01``."""
    head, peptide, *_rest = column.split("_")
    match = ALLELE_RE.match(head)
    if not match:
        return peptide, head
    gene, group, protein = match.groups()
    return peptide, f"HLA-{gene}*{group}:{protein}"


def load_features(path: Path) -> pl.DataFrame:
    frame = pl.read_csv(path)
    # Counts are floats in the dump. Binder flags are a previous call.
    drop = [c for c in frame.columns if c.endswith("_binder")]
    return frame.drop(drop)


def process_donor(dcode: Path, donor: int) -> dict[str, pl.DataFrame]:
    contig_path = dcode / f"vdj_v1_hs_aggregated_donor{donor}_all_contig_annotations.csv.gz"
    feature_path = dcode / f"vdj_v1_hs_aggregated_donor{donor}_binarized_matrix.csv.gz"
    contigs = pl.read_csv(
        contig_path,
        columns=[
            "barcode", "is_cell", "chain", "v_gene", "j_gene", "productive", "cdr3", "umis",
        ],
    )
    chains = chain_table(contigs, donor)
    edges = edge_table(chains)
    endogenous = mark_endogenous(chains)
    clones = clone_table(endogenous)
    features = load_features(feature_path)
    count_cols = [c for c in features.columns if c not in META_COLUMNS]
    dex_cols = dextramer_columns(features.columns)
    nc_cols = negative_control_columns(features.columns)
    surface = surface_columns(features.columns)
    clones = drop_doublets(clones, features, surface + dex_cols + nc_cols)
    corrected = subtract_leakage(endogenous, clones, features, count_cols)
    residuals, gammas = column_residuals(corrected, dex_cols, nc_cols)
    calls = call_specificities(clones, residuals, dex_cols, corrected, nc_cols)
    ceiling = control_envelope(clones, corrected, nc_cols)
    gamma_frame = pl.DataFrame({"column": list(gammas), "gamma": [gammas[c] for c in gammas]}).with_columns(
        pl.lit(donor).alias("donor"),
        pl.lit(ceiling).alias("envelope"),
    )
    return {
        "chains": endogenous,
        "edges": edges,
        "clones": clones,
        "calls": calls,
        "gamma": gamma_frame,
    }


def _write(frame: pl.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.write_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("--dcode", type=Path, required=True)
    parser.add_argument("--donors", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    keys = ("chains", "edges", "clones", "calls", "gamma")
    parts: dict[str, list[pl.DataFrame]] = {key: [] for key in keys}
    for donor in args.donors:
        print(f"donor {donor}", flush=True)
        result = process_donor(args.dcode, donor)
        for key, frame in result.items():
            parts[key].append(frame)
        n_clones = result["clones"]["clone_id"].n_unique() if result["clones"].height else 0
        n_ab = result["clones"].filter(pl.col("shape") == "ab").height if result["clones"].height else 0
        n_aab = result["clones"].filter(pl.col("shape") == "aab").height if result["clones"].height else 0
        print(
            f"  chains {result['chains'].height} edges {result['edges'].height} "
            f"cells ab {n_ab} aab {n_aab} clones {n_clones} calls {result['calls'].height}",
            flush=True,
        )
    sort_by = {
        "chains": ["donor", "barcode", "chain_id"],
        "edges": ["donor", "barcode", "chain_a", "chain_b"],
        "clones": ["donor", "barcode"],
        "calls": ["donor", "peptide", "mhc", "clone_id"],
        "gamma": ["donor", "column"],
    }
    for key, frames in parts.items():
        frame = pl.concat(frames, how="diagonal_relaxed")
        if frame.height:
            frame = frame.sort(sort_by[key])
        _write(frame, args.out / f"{key}.tsv")
    calls = pl.concat(parts["calls"], how="diagonal_relaxed")
    if calls.height:
        top = (
            calls.group_by("peptide", "mhc")
            .agg(pl.col("clone_id").n_unique().alias("clones"), pl.col("n_cells").sum().alias("cells"))
            .sort(["clones", "peptide", "mhc"], descending=[True, False, False])
        )
        print(top.head(15))
    else:
        print("no specificity calls")


if __name__ == "__main__":
    main()
