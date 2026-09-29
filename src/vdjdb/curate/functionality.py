"""IMGT's F / ORF / P verdict on the segment a chain names (#634).

``proofreading/imgt_alleles.tsv.gz`` is a committed, reviewed authority and three of its six columns
were read by nothing: ``functionality``, ``region_type`` and ``accession``. This module reads the
first. It answers a question no other check asks - `vdjdb qc` tests that a call *looks* like a TRBV
name and :mod:`vdjdb.curate.nomenclature` tests that IMGT *has* that name, and neither asks whether
IMGT thinks the gene is functional.

**Advisory, never a gate, and that is a statement about the biology rather than caution.** A
pseudogene V call is not automatically wrong: a P gene can rearrange, and ``TRBV21-1`` in particular
turns up in real repertoires - `annotate/junction.py` already names it as a pseudogene the
recombination model has no allele for, which is the same fact being load-bearing in one module and
invisible everywhere else. IMGT also reclassifies genes between releases, so a gate here would fail on
a reference update rather than on a curation error.

**The verdict is per species, and nothing needs deduplicating.** #634 asked for a rule to break a tie
on ``TRBV7-1*01``, which the table carries as both ``P`` and ``ORF``. Measured: that is entirely
cross-species - human ``ORF``, *Macaca fascicularis* and both *Pongo* ``P``, five other species ``F`` -
and keyed on ``(species, allele)`` the table has **zero** duplicate rows. So the resolution is a
lookup, not a reconciliation.

**The bracketed forms are kept.** IMGT writes nine: ``F``, ``(F)``, ``[F]`` and the same bracketing for
``ORF`` and ``P``, where a bracket carries whether the assignment is uncertain or allele-specific.
Collapsing them to a boolean throws that away, so the verdict ships as IMGT spells it and
:func:`is_functional` reads the letters out of it.
"""
from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path

import polars as pl

from ..config import Paths
from .nomenclature import IMGT_SPECIES

#: The three verdicts, inside whatever brackets IMGT used. ``F`` is functional; ``ORF`` is an open
#: reading frame with no evidence of expression; ``P`` is a pseudogene.
_VERDICT = re.compile(r"^[\[(]?(F|ORF|P)[\])]?$")

#: How a chain's verdict was resolved, so a reader can tell a direct answer from an inherited one.
#: ``allele`` is IMGT's verdict for the exact allele; ``gene`` is every verdict IMGT gives the alleles
#: of that gene, joined, for a call with no allele or an allele IMGT does not list; ``unlisted`` is a
#: call IMGT has at neither depth, which :mod:`vdjdb.curate.nomenclature` is the check for.
LEVELS: tuple[str, ...] = ("allele", "gene", "unlisted")

COLUMNS: tuple[str, ...] = ("record_id", "gene", "segment", "call", "functionality", "level")


@lru_cache(maxsize=4)
def _table(root: Path) -> pl.DataFrame:
    with gzip.open(root / "proofreading" / "imgt_alleles.tsv.gz") as fh:
        return pl.read_csv(fh.read(), separator="\t", infer_schema=False)


@lru_cache(maxsize=4)
def verdicts(root: Path | None = None) -> dict[tuple[str, str], tuple[str, str]]:
    """``(species, call) -> (verdict, level)``, for every allele and gene IMGT lists.

    A gene's verdict is every distinct verdict among its alleles, sorted and joined with ``/`` -
    ``F/ORF`` says the gene has functional alleles and at least one that is not. That is context for
    a reader and not a finding: see :func:`is_functional` for why ``any`` rather than ``all``.
    """
    table = _table(root or Paths.discover().root)
    out: dict[tuple[str, str], tuple[str, str]] = {}
    for vdjdb, imgt in IMGT_SPECIES.items():
        rows = table.filter(pl.col("species") == imgt)
        for allele, verdict in zip(rows["imgt_allele_id"], rows["functionality"], strict=True):
            out[(vdjdb, allele)] = (verdict, "allele")
        genes = (rows.group_by("imgt_gene_id")
                 .agg(pl.col("functionality").unique().sort().str.join("/").alias("v")))
        for gene, verdict in zip(genes["imgt_gene_id"], genes["v"], strict=True):
            out.setdefault((vdjdb, gene), (verdict, "gene"))
    return out


def is_functional(verdict: str) -> bool:
    """True where **any** verdict in ``verdict`` is a form of ``F``.

    ``any``, not ``all``, and the difference is 15,283 findings. A single-token verdict is one
    allele's, where the two readings agree. A joined one is a *gene*'s, every verdict among its
    alleles, and there the question a finding has to answer is whether the call is evidence of a
    problem - so if some allele of the gene is functional, an allele-less call on it is not.

    Measured: under ``all``, ``TRBJ2-7`` reads ``F/ORF`` because ``*02`` is an ORF, and it is one of
    the commonest J calls in VDJdb, so the advisory rule fired on **17,891 chunk rows** instead of
    2,608 chain-segments. Nothing a curator can act on is in that difference.
    """
    parts = [p for p in verdict.split("/") if p]
    return any(m.group(1) == "F" for p in parts if (m := _VERDICT.match(p)))


def resolve(species: pl.Series, calls: pl.Series,
            root: Path | None = None) -> tuple[pl.Series, pl.Series]:
    """``(verdict, level)`` per row. An empty call, or one IMGT lists nowhere, is ``("", "unlisted")``.

    Allele first, then the gene, because an allele-level verdict is the answer to the question asked
    and a gene-level one is an inference from its siblings.
    """
    table = verdicts(root)
    got = [table.get((sp, call)) or table.get((sp, call.split("*")[0])) or ("", "unlisted")
           for sp, call in zip(species, calls, strict=True)]
    return (pl.Series("functionality", [v for v, _ in got], dtype=pl.Utf8),
            pl.Series("level", [level for _, level in got], dtype=pl.Utf8))


def report(chains: pl.DataFrame, records: pl.DataFrame,
           root: Path | None = None) -> pl.DataFrame:
    """One row per chain-segment IMGT does not call functional. :data:`COLUMNS`, sorted.

    Long rather than wide - one row per (chain, segment) rather than two columns on a chain - because
    a reader groups by the call to decide what to do about it, and a wide frame makes that two
    group-bys over two columns.
    """
    keyed = chains.join(records.select("record_id", "species"), on="record_id", how="left")
    parts = []
    for segment, column in (("V", "v.segm"), ("J", "j.segm")):
        sub = keyed.filter(pl.col(column) != "")
        if sub.is_empty():
            continue
        verdict, level = resolve(sub["species"], sub[column], root)
        parts.append(sub.select("record_id", "gene", pl.lit(segment).alias("segment"),
                                pl.col(column).alias("call"))
                     .with_columns(verdict, level))
    if not parts:
        return pl.DataFrame(schema=dict.fromkeys(COLUMNS, pl.Utf8))
    return (pl.concat(parts, how="vertical")
            .filter((pl.col("level") != "unlisted")
                    & ~pl.col("functionality").map_elements(is_functional, return_dtype=pl.Boolean))
            .select(*COLUMNS)
            .sort("call", "record_id", "gene", "segment"))


def summarise(flagged: pl.DataFrame) -> pl.DataFrame:
    """``(segment, functionality, level, calls, chains)``, largest first. The advisory count."""
    if flagged.is_empty():
        return pl.DataFrame(schema={"segment": pl.Utf8, "functionality": pl.Utf8, "level": pl.Utf8,
                                    "calls": pl.Utf8, "chains": pl.UInt32})
    return (flagged.group_by("segment", "functionality", "level")
            .agg(pl.col("call").unique().sort().str.join(", ").alias("calls"),
                 pl.len().alias("chains"))
            .sort("chains", descending=True))
