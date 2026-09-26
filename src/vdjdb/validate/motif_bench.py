"""Score a motif clustering, and score legacy's the same way, on one cohort.

Phase 11's acceptance criterion is a comparison, and a comparison is only a comparison if both sides
are scored identically on identical rows. So this module does exactly two things: turn a
``cluster_members``-shaped file into the ``(label, cluster)`` frame :mod:`vdjdb.validate.metrics_lib`
expects, and restrict every clustering under comparison to one shared cohort.

⚠ **The metrics are `metrics_lib`'s, vendored verbatim, and are not interchangeable with
`mir.bench.metrics`** -- the two define `recall` over different denominators, and mixing them
silently rescales everything (ROADMAP section 8.9).

⚠ **Per-epitope clustering has purity 1.000 by construction.** A cid cannot span an epitope, so the
purity column is not evidence about a per-epitope method; it is arithmetic. Compare pooled
clusterings when purity is the question, and say which is which (ROADMAP section 8.4).

⚠ **Nothing here reads TCRvdb.** That set is held out, touched once, aggregate-only, behind
``$VDJDB_TCRVDB`` and :mod:`vdjdb.validate.guard` (ROADMAP section 11.2).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

#: The key a clustering is joined back onto -- ``cdr3|v.segm|j.segm``, the benchmark's own
#: ``production_assign`` key, and **deliberately not the epitope**. Every cid is namespaced by
#: epitope (`H.B.GILGFVFTL.1`), so joining on the epitope too would make every cluster
#: epitope-pure by construction and purity would read 1.000 for every method including legacy's --
#: measured, it does. Keying on the clonotype lets a CDR3 seen under two epitopes carry one cid,
#: which is what makes purity a number rather than arithmetic.
KEY: tuple[str, ...] = ("species", "gene", "cdr3aa", "v.segm", "j.segm")


def cohort(chains: pl.DataFrame, records: pl.DataFrame, *, species: str = "HomoSapiens",
           gene: str = "TRB", min_records: int = 30) -> pl.DataFrame:
    """One row per record in the scoring cohort: ``(species, gene, antigen.epitope, cdr3aa)``.

    Epitopes below ``min_records`` records are dropped, which is the benchmark cohort's definition.
    Every clustering compared is restricted to this, so none of them is scored on rows another one
    never saw.
    """
    df = (chains.join(records.select("record_id", "species", "antigen.epitope"), on="record_id")
                .filter((pl.col("species") == species) & (pl.col("gene") == gene)
                        & (pl.col("cdr3") != ""))
                .select("species", "gene", "antigen.epitope",
                        pl.col("cdr3").alias("cdr3aa"), "v.segm", "j.segm"))
    big = df.group_by("antigen.epitope").len().filter(pl.col("len") >= min_records).drop("len")
    return df.join(big, on="antigen.epitope").sort("antigen.epitope", *KEY)


def members_map(members: pl.DataFrame) -> pl.DataFrame:
    """``members`` reduced to one ``cluster`` per clonotype, ours or a legacy file's.

    A clonotype appearing in two clusters keeps the first in input order -- which is why the input
    order has to be the sorted one (CLAUDE.md hard rule 7). Shared with
    :func:`vdjdb.validate.qscore.frame` so the two instruments cannot disagree about which cluster a
    clonotype is in while disagreeing about how to score it.
    """
    return (members.select(*KEY, pl.col("cid").alias("cluster"))
                   .unique(subset=KEY, keep="first", maintain_order=True))


def assign(cohort_df: pl.DataFrame, members: pl.DataFrame) -> pl.DataFrame:
    """Attach ``cluster`` to every cohort record. Unclustered records get ``null``, which is what
    ``metrics_lib.precision_recall_fscore`` folds into FN.

    ``members`` is any ``cluster_members``-shaped frame -- ours or a legacy file read off disk.
    Joined on the clonotype, so a clonotype reported by several records carries its cid to all of
    them, which is how the benchmark counts.
    """
    return cohort_df.join(members_map(members), on=KEY, how="left")


def score(assigned: pl.DataFrame) -> dict:
    """``metrics_lib`` on an assigned cohort. Purity, retention, AMI, precision, recall, F1."""
    from .metrics_lib import metrics_from_assignments

    df = assigned.select("antigen.epitope", "cluster").to_pandas()
    m, _, _ = metrics_from_assignments(df)
    return m


def compare(cohort_df: pl.DataFrame, clusterings: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """Score every clustering in ``clusterings`` on ``cohort_df``. One row per method."""
    rows = []
    for name, members in clusterings.items():
        rows.append({"method": name, **score(assign(cohort_df, members))})
    return pl.DataFrame(rows)


def read_members(path: Path) -> pl.DataFrame:
    """A ``cluster_members``-shaped TSV, ours or legacy's, as an all-string frame."""
    return pl.read_csv(path, separator="\t", infer_schema_length=0, quote_char=None)


def trivial_members(cohort_df: pl.DataFrame, *, split_by_length: bool = True) -> pl.DataFrame:
    """The do-nothing partition: one cluster per epitope, every clonotype in it, nothing excluded.

    Not an algorithm -- it is what the instruments read when no clustering has happened, and it is
    the reference every candidate has to beat before its other numbers mean anything. Measured on
    the human cohort it scores ``Q`` 0.7940 / purity 0.9204 on TRA and ``Q`` 0.8632 / purity 0.9343
    on TRB, which clears three of the four admissibility axes on TRA and two on TRB
    (``docs/clustering.md`` section 8). Its lift is 1.000 by construction.

    Two things follow, and they are why this belongs in the package rather than in a sweep script:
    an absolute purity floor has to be measured against **this**, per chain per build, and ``Q`` and
    epitope coverage cannot be read as rankings.

    ``split_by_length`` is the default because **every shipped clustering is split by CDR3 length
    before it reaches a release** (``docs/clustering.md`` section 0), so the split form is the one a
    bar has to clear. It is not a detail: on human TRB the split partition scores purity 0.9343 and
    the unsplit one 0.9253, and a floor of 0.93 excludes the second while admitting the first. The
    floor the specification settled on is **0.94** on both chains, above the split form's purity on
    each (``docs/clustering.md`` section 8.1).
    """
    out = (cohort_df.select(*KEY, "antigen.epitope").unique(maintain_order=True)
           .with_columns(("trivial." + pl.col("antigen.epitope")).alias("cid")))
    if split_by_length:
        out = out.with_columns(
            (pl.col("cid") + "L" + pl.col("cdr3aa").str.len_chars().cast(pl.Utf8)).alias("cid"))
    return out


def per_epitope(cohort_df: pl.DataFrame, members: pl.DataFrame, *,
                replicated: pl.DataFrame | None = None) -> pl.DataFrame:
    """One row per epitope: how much of it is clustered, how concentrated, and whether the
    clustering tracks independent replication **there** rather than on average.

    The pooled scorecard hides the thing that matters most about a motif stage: a mean retention of
    0.33 is a different database depending on whether every epitope is a third clustered or a third
    of them are fully clustered and the rest not at all. So is a pooled lift -- it can be carried
    entirely by two large epitopes.

    Counted on **clonotypes**, not records, so a clonotype reported by forty papers does not weigh
    forty times in its own epitope's retention.

    ``percolation`` is the largest cluster's share of the epitope's clustered clonotypes: 1.0 means
    the epitope collapsed into one cluster, which is the failure mode a lift figure cannot see and
    ``Q``'s parsimony term charges for. ``singleton_clusters`` counts clusters of one, the opposite
    failure.

    Pass ``replicated`` -- ``(species, gene, cdr3aa, v.segm, j.segm, antigen.epitope, replicated)``,
    from :func:`vdjdb.motifs.tcremp.replicated` joined onto :data:`KEY` -- to get the per-epitope
    ``precision`` and ``lift``. ``lift`` is null where an epitope has no replicated clonotype at all,
    which is not a lift of zero and must not be averaged as one.

    ⚠ **Nothing here reads TCRvdb.**
    """
    clono = cohort_df.select("antigen.epitope", *KEY).unique(maintain_order=True)
    a = clono.join(members_map(members), on=KEY, how="left")

    sizes = (a.filter(pl.col("cluster").is_not_null())
              .group_by("antigen.epitope", "cluster").len().rename({"len": "csz"}))
    shape = sizes.group_by("antigen.epitope").agg(
        pl.len().alias("clusters"),
        pl.col("csz").max().alias("largest_cluster"),
        pl.col("csz").mean().round(2).alias("mean_cluster_size"),
        (pl.col("csz") == 1).sum().alias("singleton_clusters"))

    out = (a.group_by("antigen.epitope").agg(
              pl.len().alias("clonotypes"),
              pl.col("cluster").is_not_null().sum().alias("clustered"))
            .join(shape, on="antigen.epitope", how="left")
            .with_columns(pl.col("clusters", "largest_cluster", "singleton_clusters").fill_null(0))
            .with_columns(
                (pl.col("clustered") / pl.col("clonotypes")).round(4).alias("retention"),
                pl.when(pl.col("clustered") > 0)
                  .then(pl.col("largest_cluster") / pl.col("clustered"))
                  .otherwise(None).round(4).alias("percolation")))

    if replicated is not None:
        r = replicated.select("antigen.epitope", *KEY, "replicated").unique(
            subset=["antigen.epitope", *KEY], keep="first", maintain_order=True)
        j = a.join(r, on=["antigen.epitope", *KEY], how="left").with_columns(
            pl.col("replicated").fill_null(False))
        rep = j.group_by("antigen.epitope").agg(
            pl.col("replicated").sum().alias("replicated"),
            (pl.col("replicated") & pl.col("cluster").is_not_null()).sum().alias("tp"))
        out = (out.join(rep, on="antigen.epitope", how="left")
                  .with_columns(
                      pl.when(pl.col("clustered") > 0)
                        .then(pl.col("tp") / pl.col("clustered")).otherwise(None)
                        .round(4).alias("precision"))
                  .with_columns(
                      # null, not zero, where the epitope has no replicated clonotype: there is no
                      # base rate to lift against and averaging a zero in would be a fabrication.
                      pl.when((pl.col("replicated") > 0) & (pl.col("clustered") > 0))
                        .then(pl.col("precision")
                              / (pl.col("replicated") / pl.col("clonotypes")))
                        .otherwise(None).round(3).alias("lift")))
    return out.sort("clonotypes", descending=True)
