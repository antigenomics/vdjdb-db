"""The MHC dictionary: one row per distinct MHC call in the database, and its hierarchy.

Restriction is a hierarchy, and a question about it picks a level. "Is this motif HLA-A2-restricted,
A-locus-restricted, or just class I?" is three claims, so the corpus carries three MHC token families
(``m:``, ``ml:``, ``mc:``) and all three are projections of this one table. A downstream tool that
wants to group VDJdb by locus reads the mapping rather than re-deriving one, and the two cannot then
disagree.

The table is also its own proofreading report: ``status`` carries the IPD-IMGT/HLA verdict already
computed by :func:`vdjdb.assemble.epitopes.mhc_status`, so a call nothing recognises at any field depth
is visible rather than merely unmatched.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

#: Murine loci whose alleles are spelled as a haplotype suffix rather than with a ``*``: `H2-Kb` is
#: an allele of `H2-K` the way `HLA-A*02:01` is an allele of `HLA-A`, but the nomenclature gives it
#: no separator. Declared rather than pattern-matched, because the alternative is stripping a
#: trailing lowercase letter from every call and that would turn the IMGT gene `H2-Aa` into `H2-A`.
#:
#: `H2-Aa`, `H2-Ab1`, `H2-Eb1` are IMGT gene names and are their own locus, which is why they are
#: absent here (`proofreading/mhc.md`).
MURINE_LOCI: tuple[str, ...] = ("H2-IA", "H2-IE", "H2-K", "H2-D", "H2-L")

DICTIONARY_COLUMNS: tuple[str, ...] = (
    "mhc", "allele", "locus", "chain", "mhc.class", "status",
)


def two_field(call: str) -> str:
    """An allele truncated to two fields, which is the granularity VDJdb curates at.

    ``HLA-A*02:01:48`` becomes ``HLA-A*02:01``; ``HLA-A*02`` and ``H2-Kb`` are returned unchanged,
    the first because one field is all that was curated and the second because murine naming has no
    field separator. Deeper fields are not wrong, they are a resolution the database does not claim,
    and keeping them would split one restriction across several tokens.
    """
    if "*" not in call:
        return call
    gene, _, spec = call.partition("*")
    fields = spec.split(":")
    return f"{gene}*{':'.join(fields[:2])}" if fields else call


def locus(call: str) -> str:
    """The gene the call belongs to.

    Everything before the ``*`` where there is one. Otherwise the longest matching entry of
    :data:`MURINE_LOCI`, and failing that the call itself: for a murine gene name the locus and the
    molecule are the same string, because the curation does not resolve deeper. That makes ``ml:`` and
    ``m:`` the same token for such a call, which is the accurate statement rather than an invented
    hierarchy.
    """
    if "*" in call:
        return call.partition("*")[0]
    for prefix in sorted(MURINE_LOCI, key=len, reverse=True):
        if call.startswith(prefix) and call != prefix:
            return prefix
    # `B2M` and the IMGT murine gene names fall through to themselves, which is right: the class I
    # light chain is not an antigen-presenting locus and `H2-Ab1` is already a gene.
    return call


def dictionary(restriction: pl.DataFrame, *, root: Path | None = None) -> pl.DataFrame:
    """One row per distinct MHC call in ``restriction``, with its hierarchy and whether IMGT has it.

    ``chain`` is ``a`` or ``b`` as the curation placed it, so a molecule appearing on both sides of a
    class II pair gets a row per side: that placement is curated information and collapsing it would
    lose which chain a locus was reported as. ``mhc.class`` comes from the data rather than from a
    rule about the locus, because the curator's class call is the fact and a rule would be a guess
    about it.

    ``status`` is :func:`vdjdb.assemble.epitopes.mhc_status`, not a second check against the same
    IMGT table. That expression already handles the cases a fresh implementation gets wrong: a
    one-field call like ``HLA-B*07`` is a serotype-level curation and IMGT has no one-field allele
    name, so matching on full names alone would report a hundred sound calls as unknown, and murine H2
    and the light chain are outside the HLA database entirely and read ``unchecked`` rather than
    ``unknown``. Measured on the current build, exactly two calls are ``unknown``.
    """
    from ..assemble.epitopes import mhc_status

    parts = []
    for side in ("a", "b"):
        call, status = f"mhc.{side}", f"mhc.{side}.status"
        # `restriction` ships the status column, so the usual path reads it rather than recomputing
        # it. Computing it is for a frame that has not been through `assemble.epitopes`, which is
        # every test fixture.
        frame = (restriction if status in restriction.columns
                 else restriction.with_columns(mhc_status(call, root)))
        parts.append(frame.select(pl.col(call).alias("mhc"), "mhc.class",
                                 pl.lit(side).alias("chain"),
                                 pl.col(status).alias("status")))
    return (pl.concat(parts, how="vertical")
              .filter(pl.col("mhc") != "")
              .unique(subset=["mhc", "chain", "mhc.class"])
              .with_columns(
                  pl.col("mhc").map_elements(two_field, return_dtype=pl.Utf8).alias("allele"),
                  pl.col("mhc").map_elements(locus, return_dtype=pl.Utf8).alias("locus"))
              .select(list(DICTIONARY_COLUMNS))
              .sort("mhc", "chain", "mhc.class"))


def unrecognised(dict_df: pl.DataFrame) -> pl.DataFrame:
    """Calls IPD-IMGT/HLA carries at no field depth. A curation finding, one row each.

    ``unchecked`` is excluded: murine H2 and the class I light chain are outside the HLA database, so
    including them would bury two genuine findings under fifty expected ones. On the current build the
    two are ``HLA-A*08:01`` (no ``HLA-A*08`` exists at any resolution) and ``HLA-B*12`` (a serotype
    that split into B*44 and B*45 and is not a current allele group).
    """
    return dict_df.filter(pl.col("status") == "unknown").sort("mhc")
