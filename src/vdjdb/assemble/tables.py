"""The definitive tables: tidy, flat, linked by ``record_id``.

These are the pipeline's output. Everything shipped -- legacy, AIRR, the new format -- is a join and
a pivot away from here, never a parallel assembly.

One observational unit per table, one variable per column, one observation per row:

===============  ==============================  ==========================================
``records``      PK ``record_id``                 one paper's report: the pMHC, the donor and
                                                  sample, the assay, and who curated it
``chains``       PK ``(record_id, gene)``         one TCR chain
``evidence``     PK ``(record_id, evidence_id)``  one piece of supporting evidence
===============  ==============================  ==========================================

One chunk row is one record. A chunk is one publication, a row is its report on one clone, and that
row reports both chains. So ``records`` has exactly as many rows as the build reads -- 192,753 -- and
``record_id`` is a key without qualification.

``method.*`` and ``meta.*`` are on ``records`` because the README says what they are: how the
publication established the specificity. They describe the record, not the act of typing it in.
``submitter``, ``comment`` and ``chunk.id`` are the curation's own, and they sit here too, because
a record and its curation are the same row.

Identity is the complex-information columns plus the id fields -- per the README, *"duplicate
records ... will not be considered as duplicates in case they have distinct id fields"* -- plus the
chunk, because two chunks are two papers. Ids are assigned before CDR3 repair: two trimmed
sequences that repair to the same full one are still two observations, and assigning afterwards
merged 215 pairs the publications reported separately.

Three untidy shapes in the legacy format are absent here, and each of them is a join or a pivot in
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

from ..annotate.dgene import D_COLUMNS
from ..annotate.junction import NT_COLUMNS
from ..annotate.segments import INFERRED_COLUMNS
from ..identity.levels import CLONOTYPE, EPITOPE, PMHC, clone_ids, derive
from ..schema import CHAIN_COLUMNS, RECORD_COLUMNS

#: Identifies a receptor chain. Every record reporting the same chain shares one ``clonotype_id``:
#: it is the level motif evidence attaches at, and the level the independent-study support count is
#: measured on (:mod:`vdjdb.assemble.evidence`).
CLONOTYPE_KEY: tuple[str, ...] = ("species", "gene", "cdr3", "v.segm", "j.segm")

_GENES = (("alpha", "TRA"), ("beta", "TRB"))


def build_records(master: pl.DataFrame) -> pl.DataFrame:
    """One row per curated record, with ``record_id`` as a unique primary key.

    Carries the two antigen-side identifiers: ``pmhc_id`` for the peptide as presented by the
    curated alleles, ``epitope_id`` for the peptide alone. Both are derived from their own keys, so
    neither depends on what else the build contains (``ROADMAP.md`` section 10.3).
    """
    out = derive(derive(master, PMHC), EPITOPE)
    return out.select(list(RECORD_COLUMNS)).sort("record_id")


def build_chains(master: pl.DataFrame) -> pl.DataFrame:
    """One row per TCR chain: the wide alpha/beta columns pivoted long.

    A record contributes a row per chain it has, so the blanks that fill half of ``vdjdb_full.txt``
    do not exist here.
    """
    parts = []
    for gene, tag in _GENES:
        d = f"d.{gene}" if f"d.{gene}" in master.columns else None
        parts.append(
            # A chain row exists when the curator described a chain at all. A D call with no
            # CDR3 is poor curation, not absence, and dropping it would lose 34 values the legacy
            # still ships.
            master.filter((pl.col(f"cdr3.{gene}") != "")
                          | (pl.col(d) != "" if d else pl.lit(False))).select(
                pl.col("record_id"),
                pl.lit(tag).alias("gene"),
                pl.col(f"cdr3.{gene}").alias("cdr3"),
                pl.col(f"v.{gene}").alias("v.segm"),
                (pl.col(d) if d else pl.lit("")).alias("d.segm"),
                pl.col(f"j.{gene}").alias("j.segm"),
                # `species` rides along only so the clonotype key can be hashed below; it is not
                # in `CHAIN_COLUMNS`, so the final select drops it again.
                pl.col("species"),
                pl.col(f"__vend.{gene}").alias("v.end"),
                pl.col(f"__jstart.{gene}").alias("j.start"),
                pl.col(f"__cdr3old.{gene}").alias("cdr3.original"),
                pl.col(f"__fixneeded.{gene}").alias("fix.needed"),
                pl.col(f"__good.{gene}").alias("fix.good"),
                pl.col(f"__vfix.{gene}").alias("v.fix.type"),
                pl.col(f"__jfix.{gene}").alias("j.fix.type"),
                pl.col(f"__vcanon.{gene}").alias("v.canonical"),
                pl.col(f"__jcanon.{gene}").alias("j.canonical"),
                # The three stages of a segment call, side by side: what the publication reported,
                # what ships after IMGT harmonisation and allele disambiguation, and what the markup
                # engine would have called from the sequence alone. `v.segm` above is the shipped
                # one. Where the engine names no allele of its own the column is empty.
                pl.col(f"__sub.v.{gene}").alias("v.segm.submitted"),
                pl.col(f"__sub.j.{gene}").alias("j.segm.submitted"),
                (pl.col(f"__sub.{d}") if d else pl.lit("")).alias("d.segm.submitted"),
                pl.col(f"__varda.{gene}").alias("v.segm.arda"),
                pl.col(f"__jarda.{gene}").alias("j.segm.arda"),
                pl.col("TCR_hash"),
            )
        )
    # The junction-nucleotide columns are added afterwards, by `vdjdb.annotate.junction`: they need
    # the species, which is on `records`, and a model load per (species, locus).
    produced = [c for c in CHAIN_COLUMNS
                if c not in (*NT_COLUMNS, *D_COLUMNS, *INFERRED_COLUMNS, "clone_id")]
    # Sorted by the key, so the table has one order and it is the key's.
    chains = (pl.concat(parts, how="vertical")
            # 34 rows are a D call with no CDR3, so the fixer was never handed anything and left no
            # result. Empty string is the only missing marker (CLAUDE.md rule 6), `false` is what
            # "nothing was repaired" means, and -1 is what this coordinate space already reads as
            # "not mapped". The legacy export gates every one of these on `cdr3 != ""`, so nothing
            # shipped moves; a null in a shipped table is the pandas three-way ambiguity returning.
            .with_columns(
                pl.col("v.end", "j.start").fill_null(-1),
                pl.col("cdr3.original", "v.fix.type", "j.fix.type").fill_null(""),
                pl.col("fix.needed", "fix.good", "v.canonical", "j.canonical").fill_null(False),
            )
            .pipe(derive, CLONOTYPE)
            .select(produced)
            .unique(subset=["record_id", "gene"], keep="first", maintain_order=True)
            .sort("record_id", "gene"))
    # The clone keys on the pair of clonotypes, so it cannot be assigned until both exist and the
    # duplicate chain rows are gone. A record reporting one chain has no clone and gets an empty
    # string, not a null: empty is the only missing marker (hard rule 6), and a null here would be
    # the pandas three-way ambiguity coming back into a shipped table.
    return (chains.join(clone_ids(chains), on="record_id", how="left", maintain_order="left")
                  .with_columns(pl.col("clone_id").fill_null("")))


def build_tables(master: pl.DataFrame, *, release: str = "dev") -> dict[str, pl.DataFrame]:
    """The definitive tables, keyed by name."""
    from ..annotate.dgene import add_d_posterior
    from ..annotate.junction import add_junction_nt
    from ..annotate.segments import add_inferred_segments
    from ..timing import stage
    from .epitopes import build_epitopes, build_restriction
    from .evidence import build_evidence

    # Timed so the profile is an output of the build rather than something someone reconstructs with
    # an ad-hoc script afterwards. The stage names match ROADMAP_local section 49 so the numbers stay
    # comparable with the profile that produced antigenomics/vdjtools#181.
    with stage("build_records"):
        records = build_records(master)
    with stage("build_chains"):
        chains = build_chains(master)
    with stage("annotate.junction.add_junction_nt"):
        chains = add_junction_nt(chains, records)
    with stage("annotate.dgene.add_d_posterior"):
        chains = add_d_posterior(chains, records)
    with stage("annotate.segments.add_inferred_segments"):
        chains = add_inferred_segments(chains, records).select(CHAIN_COLUMNS)
    with stage("build_evidence"):
        evidence = build_evidence(records, chains, release=release)
    with stage("build_epitopes"):
        epitopes = build_epitopes(records, chains)
    with stage("build_restriction"):
        restriction = build_restriction(records)

    return {
        "records": records,
        "chains": chains,
        "evidence": evidence,
        "epitopes": epitopes,
        "restriction": restriction,
    }
