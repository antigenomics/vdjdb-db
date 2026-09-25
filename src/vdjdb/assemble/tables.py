"""The definitive tables: tidy, flat, linked by ``record_id``.

This is the pipeline's real output. Everything shipped -- legacy, AIRR, the new format -- is a
**join and a pivot** away from here, never a parallel assembly.

One observational unit per table, one variable per column, one observation per row:

===============  ==============================  ==========================================
``records``      PK ``record_id``                 one paper's report: the pMHC, the donor and
                                                  sample, the assay, and who curated it
``chains``       PK ``(record_id, gene)``         one TCR chain
``evidence``     PK ``(record_id, evidence_id)``  one piece of supporting evidence
===============  ==============================  ==========================================

**One chunk row is one record.** A chunk is one publication, a row is its report on one clone, and
that row reports both chains. So ``records`` has exactly as many rows as the build reads --
192,753 -- and ``record_id`` is a key without qualification.

``method.*`` and ``meta.*`` are on ``records`` because the README says what they are: how the
*publication* established the specificity. They describe the record, not the act of typing it in.
``submitter``, ``comment`` and ``chunk.id`` are the curation's own, and they sit here too, because
a record and its curation are the same row.

Identity is the complex-information columns plus the id fields -- per the README, *"duplicate
records ... will not be considered as duplicates in case they have distinct id fields"* -- plus the
chunk, because two chunks are two papers. Ids are assigned **before** CDR3 repair: two trimmed
sequences that repair to the same full one are still two observations, and assigning afterwards
merged 215 pairs the publications reported separately.

Three untidy shapes in the legacy format disappear here, and each of them is a join or a pivot in
the other direction:

* **paired alpha/beta columns** (``vdjdb_full.txt``) -- a chain is an observation, so it is a row.
  A record with only a beta has one row, not a row half full of blanks.
* **duplicated record fields** (``vdjdb.txt``) -- the epitope and assay are written once per record
  rather than once per chain.
* **JSON blobs** (``method``, ``meta``, ``cdr3fix``) -- every member is a column. A JSON column
  cannot be filtered, grouped or joined without parsing, and in the release the same field is a
  JSON number on one row and a string on the next.

``complex.id`` has no place here either: two chains of one clone are two rows sharing a
``record_id``. The legacy counter is assigned at export, which is the only place it means anything.
"""
from __future__ import annotations

import polars as pl

from ..schema import METHOD_COLUMNS

#: Identity and provenance, first so the tables read left to right from "which record is this".
RECORD_IDENTITY: tuple[str, ...] = ("record_id",)

#: The pMHC and its parent antigen -- what the receptor recognises.
RECORD_ANTIGEN: tuple[str, ...] = (
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species",
)

RECORD_PROVENANCE: tuple[str, ...] = ("reference.id",)

#: The donor and sample a record was observed in. Part of its identity: the same TCR against the
#: same epitope in two donors is two records, not one seen twice.
RECORD_SAMPLE: tuple[str, ...] = (
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.tissue",
)

#: Where the record was written down. Part of the record: one row, one paper, one report.
RECORD_CURATION: tuple[str, ...] = ("chunk.file", "chunk.row", "chunk.id", "submitter", "comment")

#: One row per curated record: everything the reference publication reports about it.
RECORD_COLUMNS: tuple[str, ...] = (
    *RECORD_IDENTITY, *RECORD_ANTIGEN, *RECORD_PROVENANCE, *RECORD_SAMPLE,
    "meta.epitope.id", "meta.donor.MHC", "meta.donor.MHC.method", "meta.structure.id",
    "meta.subset.frequency",
    *METHOD_COLUMNS, "method.pairing",
    "vdjdb.score",
    *RECORD_CURATION,
)

#: One row per TCR chain. ``cdr3fix.*`` is flattened: every member of the JSON blob is a column.
CHAIN_COLUMNS: tuple[str, ...] = (
    "record_id", "gene",
    "cdr3", "v.segm", "d.segm", "j.segm",
    "v.end", "j.start",
    "cdr3.original", "fix.needed", "fix.good",
    "v.fix.type", "j.fix.type", "v.canonical", "j.canonical",
    "TCR_hash",
)

_GENES = (("alpha", "TRA"), ("beta", "TRB"))


def build_records(master: pl.DataFrame) -> pl.DataFrame:
    """One row per curated record, with ``record_id`` as a unique primary key."""
    return master.select(list(RECORD_COLUMNS)).sort("record_id")


def build_chains(master: pl.DataFrame) -> pl.DataFrame:
    """One row per TCR chain: the wide alpha/beta columns pivoted long.

    A record contributes a row per chain it actually has, so the blanks that fill half of
    ``vdjdb_full.txt`` simply do not exist.
    """
    parts = []
    for gene, tag in _GENES:
        d = f"d.{gene}" if f"d.{gene}" in master.columns else None
        parts.append(
            # A chain row exists when the curator described a chain at all. A D call with no
            # CDR3 is poor curation, not absence, and dropping it would silently lose 34 values
            # the legacy still ships.
            master.filter((pl.col(f"cdr3.{gene}") != "")
                          | (pl.col(d) != "" if d else pl.lit(False))).select(
                pl.col("record_id"),
                pl.lit(tag).alias("gene"),
                pl.col(f"cdr3.{gene}").alias("cdr3"),
                pl.col(f"v.{gene}").alias("v.segm"),
                (pl.col(d) if d else pl.lit("")).alias("d.segm"),
                pl.col(f"j.{gene}").alias("j.segm"),
                pl.col(f"__vend.{gene}").alias("v.end"),
                pl.col(f"__jstart.{gene}").alias("j.start"),
                pl.col(f"__cdr3old.{gene}").alias("cdr3.original"),
                pl.col(f"__fixneeded.{gene}").alias("fix.needed"),
                pl.col(f"__good.{gene}").alias("fix.good"),
                pl.col(f"__vfix.{gene}").alias("v.fix.type"),
                pl.col(f"__jfix.{gene}").alias("j.fix.type"),
                pl.col(f"__vcanon.{gene}").alias("v.canonical"),
                pl.col(f"__jcanon.{gene}").alias("j.canonical"),
                pl.col("TCR_hash"),
            )
        )
    # Sorted by the key, so the table has one order and it is the key's.
    return (pl.concat(parts, how="vertical").select(CHAIN_COLUMNS)
            .unique(subset=["record_id", "gene"], keep="first", maintain_order=True)
            .sort("record_id", "gene"))


def build_tables(master: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """The definitive tables, keyed by name."""
    return {
        "records": build_records(master),
        "chains": build_chains(master),
    }
