"""The epitope catalogue: what VDJdb knows about each antigen, and which MHC presents it.

Two tidy tables, both derived from ``records`` and both shipped:

===============  =================================================  ===============================
``epitopes``     PK ``(antigen.epitope, antigen.species)``          the antigen itself
``restriction``  PK ``(antigen.epitope, antigen.species, mhc.a,     one presenting MHC
                 mhc.b)``
===============  =================================================  ===============================

The key is the epitope and the species, because a peptide sequence is not unique to one organism:
13 epitopes in the corpus are reported under two species, and they are not errors --
``VEALYLVCG`` is insulin B in both ``HomoSapiens``/``INS`` and ``MusMusculus``/``Ins2``, and
``KLPDDFMGC`` is conserved between SARS-CoV and SARS-CoV-2. ``patches/antigen_epitope_species_gene.dict``
is keyed on the peptide alone and cannot express any of them; this table can, which is why it is a
table rather than a view over the patch.

MHC alleles are checked, not assumed, and a name that resolves against neither authority **fails the
build**. ``proofreading/mhc_alleles.tsv.gz`` is the local mirror of IPD-IMGT/HLA
(<https://www.ebi.ac.uk/ipd/imgt/hla/>), 46,005 alleles at full four-field resolution. A VDJdb call is
usually two-field, so membership is decided by prefix: ``HLA-A*02:01`` matches the 504 rows beginning
``HLA-A*02:01:``. One-field calls such as ``HLA-A*02`` are allele groups and are accepted the same way.
``proofreading/mhc_nonhuman.tsv`` is the authority for everything IPD-IMGT/HLA cannot adjudicate --
the murine ``H2-`` names, macaque ``Mamu-A*01``, and ``B2M``, which is not an MHC gene at all.

So ``status`` reads ``known`` (IPD-IMGT/HLA), ``declared`` (the non-human table) or ``unknown``, and
``unknown`` is fatal. It used to read ``unchecked`` for a non-HLA name, meaning no authority existed;
one exists now, so the word would be false. Measured on the corpus this gate was written against:
192,793 records, **no blank MHC cell and no unresolved name in either column** -- which is what makes
failing on one a gate rather than a migration.
"""
from __future__ import annotations

import gzip
import re
from functools import lru_cache
from pathlib import Path

import polars as pl

from ..config import Paths
from ..curate.promiscuity import annotate as promiscuity_annotate
from ..schema import EPITOPE_COLUMNS

KEY: tuple[str, ...] = ("antigen.epitope", "antigen.species")

#: Not an MHC allele: the invariant light chain of every class-I molecule.
_B2M = "B2M"


#: IPD-IMGT/HLA expression suffixes, which sit on the **terminal field** of an allele name: N null,
#: L low, S secreted, Q questionable, C cytoplasmic, A aberrant.
_SUFFIX = re.compile(r"[NLSQCA]$")


@lru_cache(maxsize=4)
def _hla_prefixes(root: Path) -> frozenset[str]:
    """Every prefix of an IPD-IMGT/HLA allele name at each field depth.

    Precomputed rather than matched with a regex per call: 46,005 alleles yield a few hundred
    thousand prefixes, and membership is then a hash lookup instead of a scan.

    The terminal field is recorded both with and without its expression suffix, because the suffix is
    part of the field rather than a separate one: splitting ``HLA-A*24:09N`` on ``:`` yields
    ``HLA-A*24`` and ``HLA-A*24:09N`` and never ``HLA-A*24:09``, so a curator who writes the
    two-field name of a null allele was reported as naming nothing. 2,442 of the 53,313 prefixes are
    reachable only this way. Zero corpus strings are affected today, which is why the defect survived:
    it was waiting for the first chunk to spell one.
    """
    with gzip.open(root / "proofreading" / "mhc_alleles.tsv.gz") as fh:
        table = pl.read_csv(fh.read(), separator="\t", infer_schema=False)
    out: set[str] = set()
    for name in table["allele_name"]:
        gene, _, fields = name.partition("*")
        parts = fields.split(":")
        bare = _SUFFIX.sub("", parts[-1])
        for depth in range(1, len(parts) + 1):
            tail = parts[:depth]
            out.add(f"{gene}*{':'.join([*tail[:-1], bare if depth == len(parts) else tail[-1]])}")
        out.add(name)
    return frozenset(out)


@lru_cache(maxsize=4)
def _confirmed_prefixes(root: Path) -> frozenset[str]:
    """Every prefix under which IPD-IMGT/HLA has at least one **Confirmed** allele (#634).

    ``confirmed`` is IPD's own ``Confirmed`` / ``Unconfirmed``: whether the allele has been
    independently observed, or rests on a single submission. 13,538 of the 46,005 alleles are
    Confirmed.

    Prefix-wise and *any*, matching :func:`_hla_prefixes`, because a call is a prefix: ``HLA-A*02``
    names thousands of alleles and the question is whether any of them is confirmed, not whether all
    are.

    Measured 2026-09-29: **14 calls over 105 records**, where #634 expected four over five. It matched
    at two-field depth with the expression suffix stripped; matching the call as written finds the
    deeper ones, and those are the larger population. ``HLA-A*02:01:48`` alone is 80 records, and it
    is a third-field allele resting on one cell from one submitting group where ``HLA-A*02:01`` has 169
    Confirmed alleles under it - which is a curation question about what the submitters meant, not a
    name to reject.
    """
    with gzip.open(root / "proofreading" / "mhc_alleles.tsv.gz") as fh:
        table = pl.read_csv(fh.read(), separator="\t", infer_schema=False)
    out: set[str] = set()
    for name in table.filter(pl.col("confirmed") == "Confirmed")["allele_name"]:
        gene, _, fields = name.partition("*")
        parts = fields.split(":")
        bare = _SUFFIX.sub("", parts[-1])
        for depth in range(1, len(parts) + 1):
            tail = parts[:depth]
            out.add(f"{gene}*{':'.join([*tail[:-1], bare if depth == len(parts) else tail[-1]])}")
        out.add(name)
    return frozenset(out)


@lru_cache(maxsize=4)
def _nonhuman_names(root: Path) -> frozenset[str]:
    """The declared vocabulary for every MHC name IPD-IMGT/HLA cannot adjudicate.

    Only ``name`` is read. The table's other columns are the curator's record of what each name is;
    checking a call against the record's species would be a second gate, and this is the first.
    """
    table = pl.read_csv(root / "proofreading" / "mhc_nonhuman.tsv", separator="\t",
                        infer_schema=False, comment_prefix="#")
    return frozenset(table["name"])


def mhc_status(column: str, root: Path | None = None) -> pl.Expr:
    """``known`` / ``unconfirmed`` / ``declared`` / ``unknown`` / ``""`` for an MHC column.

    ``known`` resolves in IPD-IMGT/HLA and has at least one Confirmed allele under it;
    ``unconfirmed`` resolves there and does not; ``declared`` resolves in
    ``proofreading/mhc_nonhuman.tsv``; ``unknown`` in none of them -- which
    :func:`assert_mhc_resolves` turns into a build failure.

    **``unconfirmed`` is not fatal and is a refinement of ``known``, not a rejection** (#634). It is
    the second of two questions: #624's gate asks whether IPD carries the name at any field depth, and
    all four of these pass it. "Does anyone other than the submitter believe the allele exists" is the
    more useful question for a call in a specificity database, and an unconfirmed allele is still a
    real name - 14 calls over 105 records, 80 of them one third-field allele, so what to do about them
    is a curation decision rather than a submission error.
    """
    root = root or Paths.discover().root
    known = list(_hla_prefixes(root))
    confirmed = list(_confirmed_prefixes(root))
    declared = list(_nonhuman_names(root))
    return (
        pl.when(pl.col(column) == "").then(pl.lit(""))
        .when(pl.col(column).is_in(declared)).then(pl.lit("declared"))
        .when(pl.col(column).is_in(confirmed)).then(pl.lit("known"))
        .when(pl.col(column).is_in(known)).then(pl.lit("unconfirmed"))
        .otherwise(pl.lit("unknown"))
        .alias(f"{column}.status")
    )


def assert_mhc_resolves(records: pl.DataFrame, root: Path | None = None) -> None:
    """Fail the build on an MHC call that is blank or resolves against neither authority.

    A record whose restriction names nothing is not a record with a gap; it is a record whose
    restriction is wrong, and it has always been. Every consumer that filters VDJdb by donor type
    silently drops it, `restriction` keys on it, and `TCR_hash` hashes it, so the string reaching the
    release unchecked is worse than the build stopping. Blank is the same defect with less to go on.

    The message names the value, the column, the cell count and the chunks, because the fix is a
    `patches/mhc.dict` entry or a `proofreading/mhc_nonhuman.tsv` row and a curator needs to know
    which paper reported it.
    """
    cols = [c for c in ("mhc.a", "mhc.b") if c in records.columns]
    bad = (
        pl.concat([records.select(pl.lit(c).alias("column"), pl.col(c).alias("value"),
                                  mhc_status(c, root).alias("status"),
                                  pl.col("chunk.file") if "chunk.file" in records.columns
                                  else pl.lit("").alias("chunk.file"))
                   for c in cols], how="vertical")
        .filter(pl.col("status").is_in(["", "unknown"]))
        .group_by("column", "value")
        .agg(pl.len().alias("cells"),
             pl.col("chunk.file").unique().sort().str.join(", ").alias("chunks"))
        .sort("cells", "column", "value", descending=[True, False, False])
    )
    if bad.is_empty():
        return
    lines = [f"  {r['column']}  {r['value'] or '<blank>'!r}  {r['cells']:,} cell(s)  {r['chunks']}"
             for r in bad.iter_rows(named=True)]
    raise ValueError(
        f"{bad.height} MHC value(s) resolve against neither IPD-IMGT/HLA nor "
        f"proofreading/mhc_nonhuman.tsv:\n" + "\n".join(lines)
        + "\n\nCorrect the call in patches/mhc.dict, or declare the name in "
          "proofreading/mhc_nonhuman.tsv if it is a species the HLA database does not cover.")


#: What `proofreading/epitope_proteome.tsv` calls a peptide one residue from a proteome peptide.
#: Spelled as what was measured, not as what it might mean - see `curate.antigens`.
ONE_SUB = "one_substitution"


@lru_cache(maxsize=8)
def _analogues(root: Path) -> pl.DataFrame:
    """``antigen.epitope -> (proteome.peptide, proteome.substitution)`` from the authority.

    211 of the 791 epitopes the table covers are one substitution from a peptide in their host's
    proteome, and **67 of those have the proteome form curated in VDJdb as well** - 35,346 records,
    18.3 % of the database, sitting on rows a user cannot currently tell are related. `SLLMWITQV`
    is 29,729 records and `SLLMWITQC` is 13, so a query for NY-ESO-1 responses returns the
    anchor-optimised peptide and never learns the other row exists (#632).

    **Neither row is wrong.** The epitope sequence is the ground truth - it is the peptide the
    experiment used - and the difference from the proteome is almost always deliberate: an
    anchor-optimised vaccine peptide, a designed altered-peptide ligand, a heteroclitic variant, or
    a structure solved with a modified peptide. Of the 36,496 records on these 211 epitopes, 58
    carry a `meta.structure.id` and 54 come from `PDB_Database.tsv`, so the crystallography case is
    real and small. This column is a **link**, not a finding.

    Deduplicated on the epitope: one peptide may have a row per proteome, and `VEALYLVSG` does -
    human `INS` and mouse `Ins2` both answer `VEALYLVCG` at `8C>S`, which is the same statement
    twice rather than two.
    """
    path = root / "proofreading" / "epitope_proteome.tsv"
    schema = {"antigen.epitope": pl.String, "proteome.peptide": pl.String,
              "proteome.substitution": pl.String}
    if not path.exists():
        return pl.DataFrame(schema=schema)
    table = pl.read_csv(path, separator="\t", infer_schema=False,
                        comment_prefix="#").fill_null("")
    return (table.filter(pl.col("verdict") == ONE_SUB)
                 .select("antigen.epitope",
                         pl.col("source.peptide").alias("proteome.peptide"),
                         pl.col("source.subs").alias("proteome.substitution"))
                 .unique(subset=["antigen.epitope"], keep="first", maintain_order=True)
                 .sort("antigen.epitope"))


def build_epitopes(records: pl.DataFrame, chains: pl.DataFrame,
                   root: Path | None = None) -> pl.DataFrame:
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
            # catalogue reports the dominant one and `curate.submission.epitope_sources` counts the
            # rest into `out/reports/epitope-sources.tsv`, so the choice is visible. That comment
            # used to name an `epitopes.conflicts` that was never written, which meant the discarded
            # labels were reported nowhere -- including the 187 `Eef2`..`Eef188` labels on one
            # peptide, an index written into `antigen.gene` (#633).
            pl.col("antigen.gene").mode().sort().first().alias("antigen.gene"),
            pl.col("mhc.class").mode().sort().first().alias("mhc.class"),
            pl.len().alias("records"),
            pl.col("reference.id").n_unique().alias("references"),
        )
        .join(counts, on=KEY, how="left")
        .with_columns(pl.col("antigen.epitope").str.len_chars().cast(pl.Int64)
                      .alias("epitope.length"),
                      pl.col("chains").fill_null(0), pl.col("clonotypes").fill_null(0))
        .join(_analogues(root or Paths.discover().root), on="antigen.epitope", how="left")
        .with_columns(pl.col("proteome.peptide").fill_null(""),
                      pl.col("proteome.substitution").fill_null(""))
        .select(EPITOPE_COLUMNS)
        .sort(KEY)
    )


def build_restriction(records: pl.DataFrame, root: Path | None = None) -> pl.DataFrame:
    """One row per ``(antigen, presenting MHC)``, every call resolved or the build stops.

    The six promiscuity columns are joined on last, from the committed
    ``proofreading/epitope_promiscuity.tsv``. Nothing is predicted here: a prediction is a moving
    target, and putting one in a build would make the output depend on a model download (hard rule
    9) - `vdjdb promiscuity` refreshes that table through its own pull request.
    """
    assert_mhc_resolves(records, root)
    built = (
        records.group_by(*KEY, "mhc.a", "mhc.b", "mhc.class")
        .agg(pl.len().alias("records"),
             pl.col("reference.id").n_unique().alias("references"))
        .with_columns(mhc_status("mhc.a", root), mhc_status("mhc.b", root))
    )
    return promiscuity_annotate(built, root).sort(*KEY, "mhc.a", "mhc.b")
