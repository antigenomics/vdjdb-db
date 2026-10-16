#!/usr/bin/env python3
"""Chain co-occurrence graph for the four-donor Immudex dCODE matrices.

A node is one chain (locus + V + J + CDR3). An edge is two chains in the same
cell barcode, carrying both UMI counts. Leakage is a Poisson ambient component

    lambda(b, c) = rho_b * U_c / sum(U)

with rho_b fit once from the chains in barcode b that cannot be the primary
receptor (a third TRA, a second TRB). A chain is endogenous when that mixture
assigns it to the high component. A kept label is alpha-beta, alpha-alpha-beta,
beta only, alpha only, or two alphas. A second beta is not a label.

A barcode whose counts fit the sum of two parents is dropped. Contig UMI is not
used to rescale dextramer UMI: the multimer count on a kept barcode stays the
measured count. The clone id is only a label after that.

Three counts are not interchangeable. Cells per label are the sample size of
the within-label distribution. TCR contig UMI is the chain-graph weight only.
Dextramer UMI is the number of pMHC multimers bound: zero-inflated, median 0,
so specificity is scored on log1p(dextramer UMI).

A label is specific for a peptide when its median log1p(dextramer UMI) sits
above three distributions on that scale: the label's own per-cell negative
controls, the donor's negative-control edge (mean + 3 SD; the control median
is 0, so a median/MAD edge would be 0), and that peptide's across-label mean
+ 3 SD. A peptide that is high on most labels has a high across-label edge.
No UMI cutoff is scanned and no permutation null is used.

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
# Control median of log1p(dextramer UMI) is 0, so median and MAD mark no edge.
# Three standard deviations is the control-chart width of that distribution.
CONTROL_DEVIATIONS = 3.0
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
    """One label per barcode whose endogenous chains are a kept receptor shape.

    Kept shapes are alpha-beta, alpha-alpha-beta, beta only, alpha only, and
    two alphas. A second beta is a doublet and is not a label. The id is a
    grouping key for the dextramer distribution, not a specificity score.
    """
    empty = pl.DataFrame(
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
    keep = endogenous.filter(pl.col("endogenous"))
    if keep.is_empty():
        return empty
    tra = (
        keep.filter(pl.col("locus") == "TRA")
        .sort("chain_id")
        .group_by("donor", "barcode")
        .agg(pl.col("chain_id").alias("alphas"), pl.len().alias("n_tra"))
    )
    trb = (
        keep.filter(pl.col("locus") == "TRB")
        .sort("chain_id")
        .group_by("donor", "barcode")
        .agg(pl.col("chain_id").alias("betas"), pl.len().alias("n_trb"))
    )
    both = tra.join(trb, on=["donor", "barcode"], how="full", coalesce=True)
    both = both.with_columns(
        pl.col("n_tra").fill_null(0),
        pl.col("n_trb").fill_null(0),
    )
    both = both.filter(
        ((pl.col("n_tra") == 1) & (pl.col("n_trb") == 1))
        | ((pl.col("n_tra") == 2) & (pl.col("n_trb") == 1))
        | ((pl.col("n_tra") == 0) & (pl.col("n_trb") == 1))
        | ((pl.col("n_tra") == 1) & (pl.col("n_trb") == 0))
        | ((pl.col("n_tra") == 2) & (pl.col("n_trb") == 0))
    )
    both = both.with_columns(
        pl.when((pl.col("n_tra") == 1) & (pl.col("n_trb") == 1)).then(pl.lit("ab"))
        .when((pl.col("n_tra") == 2) & (pl.col("n_trb") == 1)).then(pl.lit("aab"))
        .when((pl.col("n_tra") == 0) & (pl.col("n_trb") == 1)).then(pl.lit("b"))
        .when((pl.col("n_tra") == 1) & (pl.col("n_trb") == 0)).then(pl.lit("a"))
        .otherwise(pl.lit("aa"))
        .alias("shape"),
        pl.col("alphas").list.get(0, null_on_oob=True).alias("alpha1"),
        pl.when(pl.col("n_tra") == 2).then(pl.col("alphas").list.get(1)).otherwise(None).alias("alpha2"),
        pl.col("betas").list.get(0, null_on_oob=True).alias("beta"),
    )
    both = both.with_columns(
        pl.when(pl.col("shape") == "ab").then(pl.col("alpha1") + pl.lit("||") + pl.col("beta"))
        .when(pl.col("shape") == "aab")
        .then(pl.col("alpha1") + pl.lit("|") + pl.col("alpha2") + pl.lit("||") + pl.col("beta"))
        .when(pl.col("shape") == "b").then(pl.lit("||") + pl.col("beta"))
        .when(pl.col("shape") == "a").then(pl.col("alpha1") + pl.lit("||"))
        .otherwise(pl.col("alpha1") + pl.lit("|") + pl.col("alpha2") + pl.lit("||"))
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


def _empty_calls() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "donor": pl.Int64,
            "clone_id": pl.String,
            "n_cells": pl.Int64,
            "column": pl.String,
            "peptide": pl.String,
            "mhc": pl.String,
            "score": pl.Float64,
        }
    )


def control_edge(nc_max: np.ndarray) -> tuple[float, float, float]:
    """Mean, SD, and mean + 3 SD of log1p(per-cell maximum negative-control UMI).

    ``nc_max`` is one dextramer-UMI number per cell, the maximum of that cell's
    negative-control reagents. It is not a TCR contig UMI and not a cell count.
    The median of this sample is 0, so the edge is the mean and SD on the
    log1p scale.
    """
    logged = np.log1p(np.asarray(nc_max, dtype=np.float64))
    if logged.size == 0:
        return 0.0, 0.0, 0.0
    mean = float(logged.mean())
    sd = float(logged.std())
    return mean, sd, mean + CONTROL_DEVIATIONS * sd


def call_specificities(
    clones: pl.DataFrame,
    counts: pl.DataFrame,
    dex_cols: list[str],
    controls: pl.DataFrame,
    nc_cols: list[str],
) -> pl.DataFrame:
    """One row per label per pMHC whose log1p dextramer median clears both distributions.

    ``counts`` are dextramer UMIs (multimers bound), after leakage subtraction.
    ``n_cells`` is stored as the within-label sample size and is not added to
    the score. The score is the median of log1p(dextramer UMI) inside the label.
    """
    if clones.is_empty() or not dex_cols:
        return _empty_calls()
    selected = counts.select(
        "barcode",
        *[pl.col(c).cast(pl.Float64).clip(lower_bound=0.0).alias(c) for c in dex_cols],
    )
    if nc_cols:
        ctrl = controls.select(
            "barcode",
            pl.max_horizontal(
                [pl.col(c).cast(pl.Float64).clip(lower_bound=0.0) for c in nc_cols]
            ).alias("_nc"),
        )
    else:
        ctrl = selected.select("barcode").unique().with_columns(pl.lit(0.0).alias("_nc"))
    joined = (
        clones.select("donor", "barcode", "clone_id")
        .join(selected, on="barcode", how="inner")
        .join(ctrl, on="barcode", how="left")
        .with_columns(pl.col("_nc").fill_null(0.0))
    )
    if joined.is_empty():
        return _empty_calls()
    _mean, _sd, fence = control_edge(joined.unique(subset=["barcode"])["_nc"].to_numpy())
    agg = joined.group_by("donor", "clone_id").agg(
        pl.len().alias("n_cells"),
        pl.col("_nc").log1p().median().alias("nc_score"),
        *[pl.col(c).log1p().median().alias(c) for c in dex_cols],
    )
    donor = agg["donor"].to_numpy()
    clone_ids = agg["clone_id"].to_list()
    n_cells = agg["n_cells"].to_numpy()
    nc_score = agg["nc_score"].to_numpy()
    rows: list[dict[str, object]] = []
    for col in dex_cols:
        peptide, mhc = parse_dextramer(col)
        score = agg[col].to_numpy().astype(np.float64)
        across = float(score.mean() + CONTROL_DEVIATIONS * score.std())
        keep = (score > nc_score) & (score > fence) & (score > across)
        for i in np.flatnonzero(keep):
            rows.append({
                "donor": int(donor[i]),
                "clone_id": clone_ids[i],
                "n_cells": int(n_cells[i]),
                "column": col,
                "peptide": peptide,
                "mhc": mhc,
                "score": float(score[i]),
            })
    if not rows:
        return _empty_calls()
    return pl.DataFrame(rows).sort("donor", "peptide", "mhc", "clone_id")


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
    dex_cols = dextramer_columns(features.columns)
    nc_cols = negative_control_columns(features.columns)
    surface = surface_columns(features.columns)
    clones = drop_doublets(clones, features, surface + dex_cols + nc_cols)
    # Dextramer UMI is the measured multimer count. Contig UMI does not rescale it.
    calls = call_specificities(clones, features, dex_cols, features, nc_cols)
    kept = features.join(clones.select("barcode").unique(), on="barcode", how="inner")
    if nc_cols and kept.height:
        nc_max = kept.select(
            pl.max_horizontal([pl.col(c).cast(pl.Float64) for c in nc_cols])
        ).to_numpy().ravel()
    else:
        nc_max = np.zeros(kept.height, dtype=np.float64)
    nc_mean, nc_sd, fence = control_edge(nc_max)
    control = pl.DataFrame({
        "donor": [donor],
        "n_cells": [int(nc_max.size)],
        "nc_mean": [nc_mean],
        "nc_sd": [nc_sd],
        "fence": [fence],
    })
    return {
        "chains": endogenous,
        "edges": edges,
        "clones": clones,
        "calls": calls,
        "control": control,
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
    keys = ("chains", "edges", "clones", "calls", "control")
    parts: dict[str, list[pl.DataFrame]] = {key: [] for key in keys}
    for donor in args.donors:
        print(f"donor {donor}", flush=True)
        result = process_donor(args.dcode, donor)
        for key, frame in result.items():
            parts[key].append(frame)
        clones = result["clones"]
        n_clones = clones["clone_id"].n_unique() if clones.height else 0
        shapes = clones["shape"].to_list() if clones.height else []
        print(
            f"  chains {result['chains'].height} edges {result['edges'].height} "
            f"cells ab {shapes.count('ab')} aab {shapes.count('aab')} "
            f"b {shapes.count('b')} a {shapes.count('a')} aa {shapes.count('aa')} "
            f"clones {n_clones} calls {result['calls'].height} "
            f"fence {float(result['control']['fence'][0]):.3f}",
            flush=True,
        )
    sort_by = {
        "chains": ["donor", "barcode", "chain_id"],
        "edges": ["donor", "barcode", "chain_a", "chain_b"],
        "clones": ["donor", "barcode"],
        "calls": ["donor", "peptide", "mhc", "clone_id"],
        "control": ["donor"],
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
