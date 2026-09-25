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


def assign(cohort_df: pl.DataFrame, members: pl.DataFrame) -> pl.DataFrame:
    """Attach ``cluster`` to every cohort record. Unclustered records get ``null``, which is what
    ``metrics_lib.precision_recall_fscore`` folds into FN.

    ``members`` is any ``cluster_members``-shaped frame -- ours or a legacy file read off disk.
    Joined on the clonotype, so a clonotype reported by several records carries its cid to all of
    them, which is how the benchmark counts.
    """
    m = (members.select(*KEY, pl.col("cid").alias("cluster"))
                .unique(subset=KEY, keep="first", maintain_order=True))
    return cohort_df.join(m, on=KEY, how="left")


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
