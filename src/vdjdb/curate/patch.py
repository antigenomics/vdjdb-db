"""Antigen nomenclature patches.

``patches/antigen_epitope_species_gene.dict`` maps an epitope to its canonical parent gene and
species. The legacy build applied it per chunk with ``chunk_df.T.apply`` -- transpose, then a
Python call per row; it is one join here.

Coverage, measured: the dictionary covers 245 epitopes, which is 154,373 of 203,308 chunk rows, and
**zero of those rows are swapped**. The ~90 records reported in #368 with ``antigen.gene`` and
``antigen.species`` transposed are therefore in the uncovered tail of 48,935 rows -- the fix is to
extend coverage, not to repair the patched set. That is phase 9.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..config import Paths

PATCH_FILE = "antigen_epitope_species_gene.dict"


def load_antigen_patch(path: Path | None = None) -> pl.DataFrame:
    """The epitope -> (gene, species) table, keyed on ``antigen.epitope``."""
    p = path or Paths.discover().patches / PATCH_FILE
    return (
        pl.read_csv(p, separator="\t", quote_char=None, infer_schema_length=0)
        .rename({"antigen.gene": "__gene", "antigen.species": "__species"})
        .select("antigen.epitope", "__gene", "__species")
        # `NA` stays a value. It is influenza *neuraminidase*, a real gene symbol, and pandas'
        # default `na_values` read it as missing -- then `if NaN:` is True in Python, so the
        # legacy assigned the missing value and blanked the gene on 14 rows.
        .fill_null("")
        # `keep="last"`, not "first": the file has eight epitopes listed twice, and pandas'
        # `.to_dict()` -- which the legacy build used -- keeps the last occurrence. Taking the
        # first instead moved antigen.gene on 1,158 rows and antigen.species on 708.
        .unique(subset="antigen.epitope", keep="last", maintain_order=True)
    )


def apply_antigen_patch(df: pl.DataFrame, patch: pl.DataFrame | None = None) -> pl.DataFrame:
    """Overwrite ``antigen.gene`` / ``antigen.species`` where the patch has a non-empty value.

    Empty patch values fall back to the curated value, which is what the legacy ``if ... else``
    did. A patch entry is an assertion about the epitope, not about the record.
    """
    p = patch if patch is not None else load_antigen_patch()
    return (
        df.join(p, on="antigen.epitope", how="left")
        .with_columns(
            pl.when(pl.col("__gene").fill_null("") != "")
            .then(pl.col("__gene")).otherwise(pl.col("antigen.gene")).alias("antigen.gene"),
            pl.when(pl.col("__species").fill_null("") != "")
            .then(pl.col("__species")).otherwise(pl.col("antigen.species")).alias("antigen.species"),
        )
        .drop("__gene", "__species")
    )
