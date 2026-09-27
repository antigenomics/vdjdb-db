"""The evidence table: one row per piece of support for a record.

Long, not wide. A record may have any number of pieces of evidence of any number of kinds -- a
motif cluster per method, a structure per PDB entry, a supporting publication per replication -- and
the wide form would be a table of mostly-empty boolean columns that grows a column per producer.
``vdjdb.parquet`` pivots it back to that shape at export (:mod:`vdjdb.emit.vdjdb3`), in the derived
direction rather than the authored one.

The first producer is independent replication, and it is the same computation as the motif tuning
objective (ROADMAP §11.1). "Two papers found this receptor against this epitope" is both the
strongest evidence the database holds about a record and the signal a clustering has to recover, so
there must not be two implementations of it that can disagree.

It is also the reason deduplication is within a chunk only. A chunk is one paper; collapsing two
chunks' matching rows into one record would delete this signal before it could be counted.

No held-out validation data is ever an evidence row. TCRvdb is proprietary, and it is read once at
the end to validate -- never to tune, never to annotate (`docs/outputs.md` §7).
"""
from __future__ import annotations

import polars as pl

from ..schema import EVIDENCE_TABLE_COLUMNS

#: The support count is per receptor chain against one epitope: ``clonotype_id`` already covers
#: species, gene, CDR3, V and J.
SUPPORT_KEY: tuple[str, ...] = ("clonotype_id", "antigen.epitope")

INDEPENDENT_STUDY = "independent_study"


def support_counts(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """Distinct reporting publications per clonotype-epitope pair, and which they are.

    One row per pair, with ``references`` the sorted distinct ``reference.id`` list and ``studies``
    its length. Phase 11 fits the clustering ``coef`` against ``studies > 1``; the evidence rows
    below are the same numbers reshaped, so the two cannot drift.

    Measured on the current corpus, human only: 4,129 of 187,238 clonotype-epitope pairs (2.21 %)
    have >= 2 distinct references when counted on `chunks/` as submitted, and 4,974 of 184,660
    when counted here, after CDR3 repair -- repair merges sequences, so pairs fall and replication
    rises. Both are the same definition at two stages of the pipeline; the tuning objective uses
    this one, because this is what ships.
    """
    return (
        chains.select("record_id", "clonotype_id")
        .join(records.select("record_id", "antigen.epitope", "reference.id"),
              on="record_id", how="left")
        .group_by(SUPPORT_KEY)
        .agg(pl.col("reference.id").unique().sort().alias("references"))
        .with_columns(pl.col("references").list.len().alias("studies"))
        .sort(SUPPORT_KEY)
    )


def independent_study(records: pl.DataFrame, chains: pl.DataFrame,
                      *, release: str = "dev") -> pl.DataFrame:
    """One evidence row per chain whose clonotype-epitope pair more than one publication reports.

    ``evidence_value`` is the other references, so a reader gets "who else saw this" rather than a
    list containing the paper they are already reading. ``evidence_score`` is the number of distinct
    references, the pair's total.
    """
    supported = support_counts(records, chains).filter(pl.col("studies") > 1)
    if supported.is_empty():
        return empty()
    return (
        chains.select("record_id", "gene", "clonotype_id")
        .join(records.select("record_id", "antigen.epitope", "reference.id"),
              on="record_id", how="left")
        .join(supported, on=SUPPORT_KEY, how="inner")
        .select(
            "record_id", "gene",
            # `(record_id, evidence_id)` is the key, so the id only has to separate a record's own
            # evidence. A readable natural value beats a hash nobody can trace back.
            pl.concat_str(pl.lit(INDEPENDENT_STUDY), pl.col("gene"), separator=":")
              .alias("evidence_id"),
            pl.lit(INDEPENDENT_STUDY).alias("evidence_type"),
            # The corpus itself is the source; there is no external database to name.
            pl.lit("").alias("evidence_source"),
            pl.col("references").list.set_difference(pl.concat_list("reference.id"))
              .list.sort().list.join(",").alias("evidence_value"),
            pl.col("studies").cast(pl.Float64).alias("evidence_score"),
            pl.lit(release).alias("first_seen_release"),
        )
        .sort("record_id", "evidence_id")
    )


def empty() -> pl.DataFrame:
    """An evidence table with no rows but the right schema, so a producer may legitimately find none."""
    return pl.DataFrame(schema={
        "record_id": pl.String, "gene": pl.String, "evidence_id": pl.String,
        "evidence_type": pl.String, "evidence_source": pl.String, "evidence_value": pl.String,
        "evidence_score": pl.Float64, "first_seen_release": pl.String,
    }).select(EVIDENCE_TABLE_COLUMNS)


def build_evidence(records: pl.DataFrame, chains: pl.DataFrame,
                   *, release: str = "dev") -> pl.DataFrame:
    """Every producer's rows, concatenated.

    One producer today. Motif clusters (phases 10-11) and structures (phase 8) append rows here
    rather than adding a column anywhere; that is what the long shape is for.
    """
    parts = [independent_study(records, chains, release=release)]
    return (pl.concat([p for p in parts if not p.is_empty()] or [empty()], how="vertical")
            .select(EVIDENCE_TABLE_COLUMNS)
            .sort("record_id", "evidence_id"))
