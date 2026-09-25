"""The epitope catalogue: what VDJdb knows about each antigen, and which MHC presents it.

Two tidy tables, both derived from ``records`` and both shipped:

===============  =================================================  ===============================
``epitopes``     PK ``(antigen.epitope, antigen.species)``          the antigen itself
``restriction``  PK ``(antigen.epitope, antigen.species, mhc.a,     one presenting MHC
                 mhc.b)``
===============  =================================================  ===============================

**Why the key is the epitope *and* the species.** A peptide sequence is not unique to one organism:
13 epitopes in the corpus are reported under two species, and they are not errors --
``VEALYLVCG`` is insulin B in both ``HomoSapiens``/``INS`` and ``MusMusculus``/``Ins2``, and
``KLPDDFMGC`` is conserved between SARS-CoV and SARS-CoV-2. ``patches/antigen_epitope_species_gene.dict``
is keyed on the peptide alone and so cannot express any of them; this table can, and that is the
reason it exists rather than being a view over the patch.

**MHC alleles are checked, not assumed.** ``proofreading/mhc_alleles.tsv.gz`` is the local mirror of
IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>), 46,005 alleles at full four-field resolution.
A VDJdb call is usually two-field, so membership is decided by **prefix**: ``HLA-A*02:01`` matches
the 504 rows beginning ``HLA-A*02:01:``. One-field calls such as ``HLA-A*02`` are allele *groups* and
are accepted the same way. Non-HLA names -- murine ``H2-Db``, ``B2M`` -- have no such authority and
are marked ``unchecked`` rather than invalid.
"""
from __future__ import annotations

import gzip
from functools import lru_cache
from pathlib import Path

import polars as pl

from ..config import Paths
from ..schema import EPITOPE_COLUMNS, RESTRICTION_COLUMNS

KEY: tuple[str, ...] = ("antigen.epitope", "antigen.species")

#: Not an MHC allele: the invariant light chain of every class-I molecule.
_B2M = "B2M"


@lru_cache(maxsize=4)
def _hla_prefixes(root: Path) -> frozenset[str]:
    """Every prefix of an IPD-IMGT/HLA allele name at each field depth.

    Precomputed rather than matched with a regex per call: 46,005 alleles yield a few hundred
    thousand prefixes, and membership is then a hash lookup instead of a scan.
    """
    with gzip.open(root / "proofreading" / "mhc_alleles.tsv.gz") as fh:
        table = pl.read_csv(fh.read(), separator="\t", infer_schema=False)
    out: set[str] = set()
    for name in table["allele_name"]:
        gene, _, fields = name.partition("*")
        parts = fields.split(":")
        for depth in range(1, len(parts) + 1):
            out.add(f"{gene}*{':'.join(parts[:depth])}")
    return frozenset(out)


def mhc_status(column: str, root: Path | None = None) -> pl.Expr:
    """``known`` / ``unknown`` / ``unchecked`` / ``""`` for an MHC column."""
    known = list(_hla_prefixes(root or Paths.discover().root))
    return (
        pl.when(pl.col(column) == "").then(pl.lit(""))
        .when(pl.col(column) == _B2M).then(pl.lit("unchecked"))
        .when(~pl.col(column).str.starts_with("HLA-")).then(pl.lit("unchecked"))
        .when(pl.col(column).str.strip_chars("N").is_in(known)).then(pl.lit("known"))
        .otherwise(pl.lit("unknown"))
        .alias(f"{column}.status")
    )


def build_epitopes(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """One row per ``(antigen.epitope, antigen.species)``, with what supports it."""
    per_record = records.select(*KEY, "antigen.gene", "mhc.class", "record_id", "reference.id")
    clono = (chains.select("record_id", "clonotype_id")
             .join(per_record.select("record_id", *KEY), on="record_id", how="inner"))
    counts = (clono.group_by(KEY)
              .agg(pl.col("clonotype_id").n_unique().alias("clonotypes"),
                   pl.len().alias("chains")))
    return (
        per_record.group_by(KEY)
        .agg(
            # One epitope under one species may still be labelled with two gene symbols; the
            # catalogue reports the dominant one and `epitopes.conflicts` in the report lists the
            # rest, rather than silently choosing.
            pl.col("antigen.gene").mode().sort().first().alias("antigen.gene"),
            pl.col("mhc.class").mode().sort().first().alias("mhc.class"),
            pl.len().alias("records"),
            pl.col("reference.id").n_unique().alias("references"),
        )
        .join(counts, on=KEY, how="left")
        .with_columns(pl.col("antigen.epitope").str.len_chars().cast(pl.Int64)
                      .alias("epitope.length"),
                      pl.col("chains").fill_null(0), pl.col("clonotypes").fill_null(0))
        .select(EPITOPE_COLUMNS)
        .sort(KEY)
    )


def build_restriction(records: pl.DataFrame, root: Path | None = None) -> pl.DataFrame:
    """One row per ``(antigen, presenting MHC)``, each allele checked against IPD-IMGT/HLA."""
    return (
        records.group_by(*KEY, "mhc.a", "mhc.b", "mhc.class")
        .agg(pl.len().alias("records"),
             pl.col("reference.id").n_unique().alias("references"))
        .with_columns(mhc_status("mhc.a", root), mhc_status("mhc.b", root))
        .select(RESTRICTION_COLUMNS)
        .sort(*KEY, "mhc.a", "mhc.b")
    )
