"""Column tuples for each table.

Phase 0 carries only what the chunk lint needs. Phase 1 (``feature/schema``) replaces this with the
full ``Field`` registry that also generates ``vdjdb.meta.txt``, the AIRR mapping and the docs tables;
until then these tuples are transcribed from ``py_src/ChunkQC.py`` and must stay identical to it.
"""
from __future__ import annotations

COMPLEX_COLUMNS: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha",
    "cdr3.beta", "v.beta", "d.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species",
    "reference.id",
)

METHOD_COLUMNS: tuple[str, ...] = (
    "method.identification", "method.frequency", "method.singlecell",
    "method.sequencing", "method.verification",
)

META_COLUMNS: tuple[str, ...] = (
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.epitope.id", "meta.tissue",
    "meta.donor.MHC", "meta.donor.MHC.method", "meta.structure.id",
)

#: Every column the build reads. `pd.concat(...)[ALL_COLS]` in the legacy pipeline drops the rest.
ALL_COLUMNS: tuple[str, ...] = COMPLEX_COLUMNS + METHOD_COLUMNS + META_COLUMNS

#: Present in real chunks, documented in the README, and silently discarded by the legacy build.
#: Listing them makes the drop deliberate and testable. ROADMAP §8 decides their fate in phase 1.
TOLERATED_DROPPED: frozenset[str] = frozenset({
    "chunk.id",
    "submitter",
    "comment",
    "meta.subset.frequency",
    "method.pairing",
})

#: Per-chunk deduplication key. Note this is *not* the score signature, which is a different
#: 11-column key that shares the name `SIGNATURE_COLS` in the legacy code. Two different keys,
#: one name -- the source of real confusion.
CHUNK_DEDUP_KEY: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha",
    "cdr3.beta", "v.beta", "d.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species", "reference.id",
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.tissue",
)

SPECIES: frozenset[str] = frozenset({
    "HomoSapiens", "MusMusculus", "RattusNorvegicus", "MacacaMulatta",
})
