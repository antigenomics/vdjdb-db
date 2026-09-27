"""The homogeneity--parsimony trade-off ``Q``, for comparing two clusterings on equal terms.

``Q = 2hp/(h+p)`` (Tiffeau-Mayer, `arXiv:2607.20799 <https://arxiv.org/abs/2607.20799>`_), the
score `~/vcs/manuscripts/2026-immrep25-audit` places its cohort ladder with. Homogeneity ``h`` asks
whether clusters predict the epitope; parsimony ``p`` asks whether they do it without shattering
each epitope into singletons. Both are normalised by the maximum attainable for the given label
partition, so a chain's epitope count and sizes divide out.

It replaces the purity/retention pair this repo used through sections 31-34, because those are two
numbers that trade against each other and ``Q`` is one that cannot be gamed from either end:

=========================  ===  ===  =====
clustering                 h    p    Q
=========================  ===  ===  =====
clusters coincide with     1    1    **1**
epitopes
every clonotype its own    1    0    **0**
cluster (shatter)
one cluster per chain      0    1    **0**
(percolate)
=========================  ===  ===  =====

Those three rows are measured against the reference implementation, and asserted in the tests.

⚠ This is the algorithm-comparison instrument, not the tuning objective. Production runs motif
detection inside one epitope's record set, so no other epitope's records are present and ``h`` is
1.000 by construction there, which is why the frame below is keyed on the clonotype, never on
``(epitope, clonotype)``. What picks the shipped configuration is the section 11.1 independent-study
objective in :mod:`vdjdb.motifs.tcremp`; see section 35.

⚠ Nothing here reads TCRvdb.
"""
from __future__ import annotations

import polars as pl

from .motif_bench import KEY, members_map


def frame(cohort_df: pl.DataFrame, members: pl.DataFrame) -> pl.DataFrame:
    """One row per ``(epitope, clonotype)`` with the cluster it landed in.

    A clonotype recorded against two epitopes appears twice with the same cluster both times,
    which is what lets ``h`` fall below 1. Unclustered clonotypes are kept as singletons rather than
    dropped: that is how ``p`` charges a clustering for low coverage, and why retention does not
    need to be a separate axis.
    """
    return (cohort_df.select("antigen.epitope", *KEY).unique(maintain_order=True)
            .join(members_map(members), on=KEY, how="left")
            .with_row_index("__i")
            .with_columns(pl.col("cluster").fill_null("__singleton." + pl.col("__i").cast(pl.Utf8)))
            .drop("__i"))


def score(cohort_df: pl.DataFrame, members: pl.DataFrame) -> dict:
    """``Q``, ``h`` and ``p`` for one clustering on one cohort, with the counts behind them."""
    import clustereval

    f = frame(cohort_df, members)
    y, k = f["antigen.epitope"].to_numpy(), f["cluster"].to_numpy()
    h = float(clustereval.homogeneity_score(y, k))
    p = float(clustereval.parsimony_score(y, k))
    real = f.filter(~pl.col("cluster").str.starts_with("__singleton."))
    return {"q": 0.0 if h + p == 0 else 2 * h * p / (h + p), "h": h, "p": p,
            "pairs": f.height, "clustered": real.height, "clusters": real["cluster"].n_unique()}


def compare(cohort_df: pl.DataFrame, clusterings: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """One row per clustering, ranked by ``Q``."""
    return (pl.DataFrame([{"method": n, **score(cohort_df, m)} for n, m in clusterings.items()])
            .sort("q", descending=True))
