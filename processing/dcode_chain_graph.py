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
the within-label likelihood. TCR contig UMI is the chain-graph weight only.
Dextramer UMI is the number of pMHC multimers bound.

Specificity is a per-peptide negative-binomial mixture. A cell's mean is a
peptide baseline times a capture scale. The scale is the posterior mean of a
gamma factor given that cell's negative-control total, and it is 1 when the
controls are zero or have no extra-Poisson variance. Background cells share
one baseline. A specific label draws one shared baseline from a discrete
uniform prior on the log baseline, supported only above that peptide's
background baseline. A label with two or more cells can instead be one
elevated cell and the rest background. Every hypothesis uses that dispersion.
Given that a label is not background, the prior splits equally between the
shared mean and the single hot cell. The call is the maximum a posteriori
hypothesis, and a tie stays background.
Median and MAD of log1p(dextramer UMI) are written afterwards as checks. They
are not inputs.

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
# Dispersion floor keeps a degenerate infinite-variance negative binomial from
# making every count look like background. It is not a count cutoff.
THETA_FLOOR = 0.05
THETA_CEILING = 1.0e4
MU_FLOOR = 1.0e-4
GRID_SIZE = 96
EM_ITERATIONS = 40
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


def _empty_checks() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "donor": pl.Int64,
            "peptide": pl.String,
            "mhc": pl.String,
            "mu0": pl.Float64,
            "theta": pl.Float64,
            "pi_background": pl.Float64,
            "pi_specific": pl.Float64,
            "pi_outlier": pl.Float64,
            "n_labels": pl.Int64,
            "n_called": pl.Int64,
            "median_log1p_called": pl.Float64,
            "mad_log1p_called": pl.Float64,
            "median_log1p_rest": pl.Float64,
            "mad_log1p_rest": pl.Float64,
            "median_log1p_nc_called": pl.Float64,
            "mad_log1p_nc_called": pl.Float64,
            "median_log1p_nc_rest": pl.Float64,
            "mad_log1p_nc_rest": pl.Float64,
        }
    )


def _logsumexp(values: np.ndarray, axis: int = 0) -> np.ndarray:
    peak = np.max(values, axis=axis, keepdims=True)
    safe = np.where(np.isfinite(peak), peak, 0.0)
    total = np.log(np.sum(np.exp(values - safe), axis=axis))
    return total + np.squeeze(safe, axis=axis)


def _segment_logsumexp(values: np.ndarray, starts: np.ndarray) -> np.ndarray:
    peak = np.maximum.reduceat(values, starts)
    safe = np.where(np.isfinite(peak), peak, 0.0)
    lengths = np.diff(np.append(starts, values.size))
    shifted = np.exp(values - np.repeat(safe, lengths))
    return np.log(np.add.reduceat(shifted, starts)) + safe


def _clone_kernel(
    n: np.ndarray,
    total: np.ndarray,
    mean: np.ndarray | float,
    theta: float,
) -> np.ndarray:
    """Mean-dependent negative-binomial log likelihood of a shared mean.

    ``n`` and ``total`` are the cell count and the summed dextramer UMI of each
    label. Gamma terms do not depend on the mean, so they cancel in every
    Bayes factor and are left out.
    """
    log_theta = math.log(theta)
    mean_arr = np.asarray(mean, dtype=np.float64)
    if mean_arr.ndim == 0:
        loc = float(mean_arr)
        return (
            n * theta * log_theta
            - (n * theta + total) * math.log(theta + loc)
            + np.where(total > 0.0, total * math.log(loc), 0.0)
        )
    return (
        n[:, None] * theta * log_theta
        - (n[:, None] * theta + total[:, None]) * np.log(theta + mean_arr)[None, :]
        + np.where(total[:, None] > 0.0, total[:, None] * np.log(mean_arr)[None, :], 0.0)
    )


def _elevated_log_prior(mu_grid: np.ndarray, mu0: float) -> np.ndarray:
    """Discrete uniform prior on log mean, supported only above the background mean."""
    log_prior = np.zeros(mu_grid.shape, dtype=np.float64)
    log_prior[mu_grid <= mu0] = -np.inf
    if not np.isfinite(log_prior).any():
        return log_prior
    return log_prior - _logsumexp(log_prior)


def _log1p_median_mad(values: np.ndarray) -> tuple[float | None, float | None]:
    """Median and MAD of log1p(count). A check, never an assignment rule."""
    if values.size == 0:
        return None, None
    logged = np.log1p(np.clip(values, 0.0, None))
    median = float(np.median(logged))
    mad = float(np.median(np.abs(logged - median)))
    return median, mad


def _stickiness(nc: np.ndarray) -> np.ndarray:
    """Per-cell scale for the peptide mean, from that cell's negative-control total.

    The prior on the scale has weight at least one observation at the control
    mean, so one extreme control cell cannot drive every other scale to zero.
    All-zero controls carry no stickiness and return scale 1.
    """
    if nc.size == 0 or float(np.sum(nc)) <= 0.0:
        return np.ones(nc.size, dtype=np.float64)
    mu = float(nc.mean())
    var = float(nc.var())
    if var <= mu + 1e-8:
        return np.ones(nc.size, dtype=np.float64)
    alpha = max(mu * mu / (var - mu), mu)
    return (alpha + nc) / (alpha + mu)


def _nb_mean_offset(y: np.ndarray, sigma: np.ndarray, theta: float, mu_init: float) -> float:
    """Negative-binomial MLE of a baseline mean with per-cell scales."""
    mu = max(mu_init, MU_FLOOR)
    for _ in range(6):
        mean = mu * sigma
        slope = np.sum(y / mu - (theta + y) * sigma / (theta + mean))
        curve = np.sum(-y / mu**2 + (theta + y) * (sigma**2) / (theta + mean) ** 2)
        if not np.isfinite(slope) or not np.isfinite(curve) or curve >= -1e-12:
            break
        nxt = mu - slope / curve
        if nxt <= 0.0:
            nxt = mu * 0.5
        if abs(nxt - mu) <= 1e-4 * max(mu, MU_FLOOR):
            return max(float(nxt), MU_FLOOR)
        mu = float(nxt)
    return max(float(mu), MU_FLOOR)


def _fit_peptide_offset(
    y: np.ndarray,
    starts: np.ndarray,
    lengths: np.ndarray,
    sigma: np.ndarray,
) -> dict[str, object]:
    """Same three hypotheses. Each cell mean is the baseline times its scale.

    The likelihood of a cell depends on that cell only through its count and
    its scale, so identical pairs are evaluated once.
    """
    n_labels = int(lengths.size)
    n = lengths.astype(np.float64)
    sigma = np.maximum(sigma, MU_FLOOR)
    uniq, inverse = np.unique(np.column_stack((y, sigma)), axis=0, return_inverse=True)
    uy = uniq[:, 0]
    us = uniq[:, 1]
    hi = math.log(max(float(np.max(uy / us)) * 1.5, 1.0))
    mu_grid = np.exp(np.linspace(-6.0, hi, GRID_SIZE))
    mu0 = max(float(np.sum(y) / np.sum(sigma)), MU_FLOOR)
    theta = 5.0
    pi = np.array([0.94, 0.03, 0.03])

    def e_step(mu0: float, theta: float, pi: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        log_prior = _elevated_log_prior(mu_grid, mu0)
        log_theta = math.log(theta)
        mean = us[:, None] * mu_grid[None, :]
        ll_u = (
            theta * log_theta
            - (theta + uy[:, None]) * np.log(theta + mean)
            + np.where(uy[:, None] > 0.0, uy[:, None] * np.log(mean), 0.0)
        )
        ll_grid = ll_u[inverse]
        clone_ll = np.add.reduceat(ll_grid, starts, axis=0)
        log_m1 = _logsumexp(ll_u + log_prior[None, :], axis=1)[inverse]
        mean0 = us * mu0
        ll0 = (
            theta * log_theta
            - (theta + uy) * np.log(theta + mean0)
            + np.where(uy > 0.0, uy * np.log(mean0), 0.0)
        )[inverse]
        clone_ll0 = np.add.reduceat(ll0, starts)
        log_m_spec = _logsumexp(clone_ll + log_prior[None, :], axis=1)
        log_m_out = _segment_logsumexp(log_m1 - ll0, starts) + clone_ll0 - np.log(n)
        log_post = np.column_stack((
            clone_ll0 + math.log(pi[0]),
            log_m_spec + math.log(pi[1]),
            log_m_out + math.log(pi[2]),
        ))
        # A singleton cannot be one hot cell plus background, so the outlier
        # prior mass belongs to the specific hypothesis.
        singleton = lengths == 1
        log_post[singleton, 1] = log_m_spec[singleton] + math.log(pi[1] + pi[2])
        log_post[singleton, 2] = -np.inf
        post = np.exp(log_post - _logsumexp(log_post, axis=1)[:, None])
        return post, log_m1, ll0

    for iteration in range(EM_ITERATIONS):
        post, log_m1, ll0 = e_step(mu0, theta, pi)
        excess = log_m1 - ll0
        peak = np.maximum.reduceat(excess, starts)
        safe = np.where(np.isfinite(peak), peak, 0.0)
        shifted = np.exp(excess - np.repeat(safe, lengths))
        segment_sum = np.add.reduceat(shifted, starts)
        p_cell = shifted / np.repeat(np.maximum(segment_sum, 1e-300), lengths)
        weight = np.repeat(post[:, 0], lengths) + np.repeat(post[:, 2], lengths) * (1.0 - p_cell)
        weight = np.where(np.isfinite(weight), weight, 0.0)
        background = weight >= 0.5
        mu_new = mu0
        if int(background.sum()) >= 2:
            y_bg = y[background]
            sigma_bg = sigma[background]
            guess = max(float(np.sum(y_bg) / np.sum(sigma_bg)), MU_FLOOR)
            mu_new = _nb_mean_offset(y_bg, sigma_bg, theta, guess)
            fitted = mu_new * sigma_bg
            numer = float(np.sum(fitted * fitted))
            denom = float(np.sum((y_bg - fitted) ** 2 - fitted))
            if denom > 1e-8 and numer > 0.0:
                theta = float(np.clip(numer / denom, THETA_FLOOR, THETA_CEILING))
            else:
                theta = THETA_CEILING
        pi_new = _hypothesis_weights(post, n_labels)
        stable = (
            iteration >= 5
            and abs(mu_new - mu0) <= 1e-3 * max(mu0, MU_FLOOR)
            and float(np.max(np.abs(pi_new - pi))) < 1e-4
        )
        mu0 = mu_new
        pi = pi_new
        if stable:
            break
    post, _, _ = e_step(mu0, theta, pi)
    call = (post[:, 1] > post[:, 0]) & (post[:, 1] > post[:, 2])
    return {"call": call, "posterior": post[:, 1], "mu0": float(mu0), "theta": float(theta), "pi": pi}


def _fit_peptide(
    y: np.ndarray,
    starts: np.ndarray,
    lengths: np.ndarray,
    sigma: np.ndarray | None = None,
) -> dict[str, object]:
    """MAP labels for one peptide. ``y`` is dextramer UMI grouped by label.

    A constant scale, including the all-ones scale of a zero control column,
    uses the shared-mean update. A varying scale is a per-cell multiplier.
    """
    n_labels = int(lengths.size)
    none = {
        "call": np.zeros(n_labels, dtype=bool),
        "posterior": np.zeros(n_labels, dtype=np.float64),
        "mu0": 0.0,
        "theta": None,
        "pi": np.array([1.0, 0.0, 0.0]),
    }
    if y.size == 0 or float(y.max()) <= 0.0:
        return none
    if (
        sigma is not None
        and sigma.size == y.size
        and float(np.max(sigma) - np.min(sigma)) > 1e-8
    ):
        return _fit_peptide_offset(y, starts, lengths, sigma)
    n = lengths.astype(np.float64)
    total = np.add.reduceat(y, starts)
    hi = math.log(max(float(y.max()) * 1.5, 1.0))
    log_mu = np.linspace(-6.0, hi, GRID_SIZE)
    mu_grid = np.exp(log_mu)
    mu0 = max(float(y.mean()), MU_FLOOR)
    theta = 5.0
    pi = np.array([0.94, 0.03, 0.03])

    def e_step(mu0: float, theta: float, pi: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        log_prior = _elevated_log_prior(mu_grid, mu0)
        clone_ll = _clone_kernel(n, total, mu_grid, theta)
        unique_y, inverse = np.unique(y, return_inverse=True)
        cell_ll = _clone_kernel(np.ones(unique_y.size), unique_y, mu_grid, theta)
        log_m1 = _logsumexp(cell_ll + log_prior[None, :], axis=1)[inverse]
        ll0 = _clone_kernel(np.ones(y.size), y, mu0, theta)
        clone_ll0 = np.add.reduceat(ll0, starts)
        log_m_spec = _logsumexp(clone_ll + log_prior[None, :], axis=1)
        log_m_out = _segment_logsumexp(log_m1 - ll0, starts) + clone_ll0 - np.log(n)
        log_post = np.column_stack((
            clone_ll0 + math.log(pi[0]),
            log_m_spec + math.log(pi[1]),
            log_m_out + math.log(pi[2]),
        ))
        # A singleton cannot be one hot cell plus background, so the outlier
        # prior mass belongs to the specific hypothesis.
        singleton = lengths == 1
        log_post[singleton, 1] = log_m_spec[singleton] + math.log(pi[1] + pi[2])
        log_post[singleton, 2] = -np.inf
        post = np.exp(log_post - _logsumexp(log_post, axis=1)[:, None])
        return post, log_m1, ll0

    for iteration in range(EM_ITERATIONS):
        post, log_m1, ll0 = e_step(mu0, theta, pi)
        excess = log_m1 - ll0
        peak = np.maximum.reduceat(excess, starts)
        safe = np.where(np.isfinite(peak), peak, 0.0)
        shifted = np.exp(excess - np.repeat(safe, lengths))
        segment_sum = np.add.reduceat(shifted, starts)
        p_cell = shifted / np.repeat(np.maximum(segment_sum, 1e-300), lengths)
        weight = np.repeat(post[:, 0], lengths) + np.repeat(post[:, 2], lengths) * (1.0 - p_cell)
        weight = np.where(np.isfinite(weight), weight, 0.0)
        # Soft mean, so one elevated cell does not set the background level.
        # Dispersion uses only cells the model explains by the background, and
        # that set's own mean: a soft mean under a hard set makes a single 1
        # look infinitely overdispersed and then no count is evidence.
        weight_sum = float(weight.sum())
        mu_new = mu0
        if weight_sum > 0.0:
            mu_new = max(float(np.sum(weight * y) / weight_sum), MU_FLOOR)
        background = weight >= 0.5
        if int(background.sum()) >= 2:
            y_bg = y[background]
            mu_bg = float(y_bg.mean())
            var = float(np.mean((y_bg - mu_bg) ** 2))
            if var > mu_bg + 1e-8 and mu_bg > 0.0:
                theta = float(np.clip(mu_bg ** 2 / (var - mu_bg), THETA_FLOOR, THETA_CEILING))
            else:
                theta = THETA_CEILING
        pi_new = _hypothesis_weights(post, n_labels)
        stable = (
            iteration >= 5
            and abs(mu_new - mu0) <= 1e-3 * max(mu0, MU_FLOOR)
            and float(np.max(np.abs(pi_new - pi))) < 1e-4
        )
        mu0 = mu_new
        pi = pi_new
        if stable:
            break
    post, _, _ = e_step(mu0, theta, pi)
    call = (post[:, 1] > post[:, 0]) & (post[:, 1] > post[:, 2])
    return {"call": call, "posterior": post[:, 1], "mu0": mu0, "theta": theta, "pi": pi}


def _hypothesis_weights(post: np.ndarray, n_labels: int) -> np.ndarray:
    """Background weight is empirical Bayes. The other two weights stay equal.

    Before the counts are seen, a one-hot label and a coherent label are
    equally likely. A collapsed outlier weight turns the one-hot label into
    a specific call.
    """
    pi_background = (float(post[:, 0].sum()) + 0.5) / (n_labels + 1.5)
    rest = (1.0 - pi_background) / 2.0
    return np.array([pi_background, rest, rest])


def _label_segments(clone_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if clone_id.size == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    change = np.flatnonzero(clone_id[1:] != clone_id[:-1]) + 1
    starts = np.concatenate((np.zeros(1, dtype=np.int64), change.astype(np.int64)))
    lengths = np.diff(np.append(starts, clone_id.size))
    return starts, lengths


def assign_specificities(
    clones: pl.DataFrame,
    counts: pl.DataFrame,
    dex_cols: list[str],
    controls: pl.DataFrame,
    nc_cols: list[str],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Assign each label by the posterior mode, and record median/MAD checks.

    ``counts`` are raw dextramer UMIs. Contig UMI is not a column this reads.
    The negative-control sum sets each cell's capture scale. The negative-control
    maximum is only the check. ``score`` on a called row is the posterior
    probability of the specific hypothesis.
    """
    if clones.is_empty() or not dex_cols:
        return _empty_calls(), _empty_checks()
    selected = counts.select(
        "barcode",
        *[pl.col(c).cast(pl.Float64).clip(lower_bound=0.0).alias(c) for c in dex_cols],
    )
    if nc_cols:
        nc_terms = [pl.col(c).cast(pl.Float64).clip(lower_bound=0.0) for c in nc_cols]
        ctrl = controls.select(
            "barcode",
            pl.max_horizontal(nc_terms).alias("_nc"),
            pl.sum_horizontal(nc_terms).alias("_nc_sum"),
        )
    else:
        ctrl = selected.select("barcode").unique().with_columns(
            pl.lit(0.0).alias("_nc"),
            pl.lit(0.0).alias("_nc_sum"),
        )
    joined = (
        clones.select("donor", "barcode", "clone_id")
        .join(selected, on="barcode", how="inner")
        .join(ctrl, on="barcode", how="left")
        .with_columns(pl.col("_nc").fill_null(0.0))
    )
    if joined.is_empty():
        return _empty_calls(), _empty_checks()
    rows: list[dict[str, object]] = []
    checks: list[dict[str, object]] = []
    for donor in joined["donor"].unique().sort().to_list():
        part = joined.filter(pl.col("donor") == donor).sort("clone_id")
        clone_id = part["clone_id"].to_numpy()
        starts, lengths = _label_segments(clone_id)
        nc = part["_nc"].to_numpy().astype(np.float64)
        sigma = _stickiness(part["_nc_sum"].to_numpy().astype(np.float64))
        label_id = clone_id[starts]
        label_n = lengths.astype(np.int64)
        for col in dex_cols:
            peptide, mhc = parse_dextramer(col)
            y = part[col].to_numpy().astype(np.float64)
            fit = _fit_peptide(y, starts, lengths, sigma)
            call = fit["call"]
            posterior = fit["posterior"]
            called_cells = np.repeat(call, lengths)
            med_c, mad_c = _log1p_median_mad(y[called_cells])
            med_r, mad_r = _log1p_median_mad(y[~called_cells])
            nc_c, nc_mad_c = _log1p_median_mad(nc[called_cells])
            nc_r, nc_mad_r = _log1p_median_mad(nc[~called_cells])
            pi = fit["pi"]
            checks.append({
                "donor": int(donor),
                "peptide": peptide,
                "mhc": mhc,
                "mu0": float(fit["mu0"]),
                "theta": None if fit["theta"] is None else float(fit["theta"]),
                "pi_background": float(pi[0]),
                "pi_specific": float(pi[1]),
                "pi_outlier": float(pi[2]),
                "n_labels": int(lengths.size),
                "n_called": int(np.sum(call)),
                "median_log1p_called": med_c,
                "mad_log1p_called": mad_c,
                "median_log1p_rest": med_r,
                "mad_log1p_rest": mad_r,
                "median_log1p_nc_called": nc_c,
                "mad_log1p_nc_called": nc_mad_c,
                "median_log1p_nc_rest": nc_r,
                "mad_log1p_nc_rest": nc_mad_r,
            })
            for i in np.flatnonzero(call):
                rows.append({
                    "donor": int(donor),
                    "clone_id": label_id[i],
                    "n_cells": int(label_n[i]),
                    "column": col,
                    "peptide": peptide,
                    "mhc": mhc,
                    "score": float(posterior[i]),
                })
    call_frame = pl.DataFrame(rows).sort("donor", "peptide", "mhc", "clone_id") if rows else _empty_calls()
    check_frame = pl.DataFrame(checks).sort("donor", "peptide", "mhc") if checks else _empty_checks()
    return call_frame, check_frame


def call_specificities(
    clones: pl.DataFrame,
    counts: pl.DataFrame,
    dex_cols: list[str],
    controls: pl.DataFrame,
    nc_cols: list[str],
) -> pl.DataFrame:
    """Called labels only. ``score`` is the posterior of the specific hypothesis."""
    calls, _checks = assign_specificities(clones, counts, dex_cols, controls, nc_cols)
    return calls


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
    calls, checks = assign_specificities(clones, features, dex_cols, features, nc_cols)
    return {
        "chains": endogenous,
        "edges": edges,
        "clones": clones,
        "calls": calls,
        "checks": checks,
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
    keys = ("chains", "edges", "clones", "calls", "checks")
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
            f"clones {n_clones} calls {result['calls'].height}",
            flush=True,
        )
    sort_by = {
        "chains": ["donor", "barcode", "chain_id"],
        "edges": ["donor", "barcode", "chain_a", "chain_b"],
        "clones": ["donor", "barcode"],
        "calls": ["donor", "peptide", "mhc", "clone_id"],
        "checks": ["donor", "peptide", "mhc"],
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
