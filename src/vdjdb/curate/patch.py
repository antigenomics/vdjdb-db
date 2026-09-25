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

#: Epitopes the patch lists **twice with different answers**. ``keep="last"`` then makes file order
#: decide, silently, which is why they are named here: the set is asserted in the tests, so it can
#: shrink but not grow.
#:
#: Two are pure nomenclature -- the EBNA3 family carries a modern and a legacy name for the same
#: protein (EBNA3B = EBNA4, EBNA3C = EBNA6) -- and the rest are genuine disagreements a curator has
#: to settle. The three HIV-1 ones sit in the Gag-Pol frameshift region, where both readings are
#: defensible and the right answer depends on which frame the study assayed.
CONFLICTING_EPITOPES: dict[str, tuple[str, ...]] = {
    "VTEHDTLLY": ("CMV IE1", "CMV pp50"),
    "EENLLDFVRF": ("EBV EBNA3A", "EBV EBNA6"),
    "IVTDFSVIK": ("EBV EBNA3B", "EBV EBNA4"),        # the same protein under two names
    "TAFTIPSI": ("HIV-1 Gag", "HIV-1 Pol"),          # Gag-Pol frameshift
    "CTPYDINQM": ("HIV-1 Gag", "HIV-1 Pol"),
    "ISPRTLNAW": ("HIV-1 Gag", "HIV-1 Pol"),
}


def conflicts(path: Path | None = None) -> pl.DataFrame:
    """Epitopes the patch answers more than one way. Empty is the goal; the build reports it."""
    p = path or Paths.discover().patches / PATCH_FILE
    raw = pl.read_csv(p, separator="\t", quote_char=None, infer_schema_length=0).fill_null("")
    return (raw.group_by("antigen.epitope")
            .agg(pl.concat_str("antigen.species", "antigen.gene", separator=" ")
                 .unique().sort().alias("answers"))
            .filter(pl.col("answers").list.len() > 1)
            .sort("antigen.epitope"))


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


def render_patch_renames(reference: pl.DataFrame, patched: pl.DataFrame) -> str:
    """``[[rename]]`` lines for the antigen patch, scoped to ``vdjdb.slim.txt``.

    ``antigen.gene`` and ``antigen.species`` are part of the slim **grouping** key, so correcting a
    gene symbol splits and merges slim rows -- 10,563 records' worth. The other two files carry the
    correction as a plain changed cell and need no rename.

    ``reference`` is the released table, **not** the raw chunks: the release was built with the patch
    as it stood then, so the only differences that reach the ledger are the entries added since. A
    rename generated against the chunks would declare hundreds of corrections the release already
    carries, and every one of them would be stale.

    The predicate is ``when_equals`` on the epitope, not ``when_contains``: a 9-mer epitope really is
    a substring of a 10-mer one (``SPRWYFYYL`` inside ``LSPRWYFYYL``), so a substring test would
    rewrite the wrong rows.
    """
    import json

    def answers(df: pl.DataFrame) -> pl.DataFrame:
        return df.select("antigen.epitope", "antigen.gene", "antigen.species").unique()

    # One answer is required on the **candidate** side only. The reference is allowed several: an
    # epitope labelled both `S` and `Spike` is exactly what the patch merges, and a rename from each
    # reference value onto the single candidate one is what makes slim's grouping line up.
    target = (answers(patched).group_by("antigen.epitope")
              .agg(pl.all().sort_by("antigen.gene")))
    target = (target.filter(pl.col("antigen.gene").list.len() == 1)
              .select("antigen.epitope",
                      pl.col("antigen.gene").list.first().alias("gene1"),
                      pl.col("antigen.species").list.first().alias("species1")))
    joined = answers(reference).join(target, on="antigen.epitope", how="inner")

    out = []
    for column, b, a in (("antigen.gene", "antigen.gene", "gene1"),
                         ("antigen.species", "antigen.species", "species1")):
        changed = (joined.filter(pl.col(b) != pl.col(a))
                   .select("antigen.epitope", b, a).unique().sort("antigen.epitope", b))
        for epitope, old_v, new_v in changed.iter_rows():
            out += ["[[rename]]",
                    f"columns = {json.dumps(column)}",
                    f"from = {json.dumps(old_v)}",
                    f"to = {json.dumps(new_v)}",
                    'when_columns = "antigen.epitope"',
                    f"when_equals = {json.dumps(epitope)}",
                    'files = "vdjdb.slim.txt"',
                    ""]
    return "\n".join(out)
