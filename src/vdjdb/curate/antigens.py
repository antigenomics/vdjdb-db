"""Where each self epitope occurs in its own species' reference proteome, and how exactly. #632.

**What this is and is not.** The antigen a receptor engages is the **pMHC**, and the epitope is the
part of it that is recognised (``docs/standards/terminology.md``). ``antigen.gene`` and
``antigen.species`` are neither: they are the peptide's **provenance**, which is why they are
deliberately absent from :data:`vdjdb.compare.diff.KEYS` while `restriction` is keyed on
``(epitope, species, mhc.a, mhc.b)`` and ``pmhc_id`` identifies the complex.

So nothing here bears on what any receptor recognises, and a disagreement in this table is never a
correctness failure in a specificity claim. What provenance *is* good for is querying - "which
receptors were found against human neoantigens" is exactly the question these columns answer - and a
label that is a protein name where a gene symbol belongs makes that query quietly incomplete. That is
the whole of what this module is for.

Within that bound, it is the one curated field with no authority at all. V/D/J calls are checked
against `imgt_alleles.tsv.gz`, MHC calls against `mhc_alleles.tsv.gz` and `mhc_nonhuman.tsv` -
fatally, since #624 - and references against `reference_ids.tsv`. ``antigen.gene`` is 427 distinct
free-text values over 2,131 epitopes, checked against a hand-kept alias list and nothing else.

The authority is already a core dependency. `mhcmatch.proteome.Proteome` carries the UniProt
reference proteomes, and this asks it one question per epitope: is this peptide in that proteome, and
if not, what is the nearest peptide that is?

**Three verdicts, and only one of them is a finding.**

`exact`
    The peptide is in the proteome at a named protein and position, and the proteome's own ``GN=``
    field gives the gene symbol - the authority `antigen.gene` has never had. 315 of 746
    curated-human epitopes and 15 of 44 curated-mouse ones. Where that symbol disagrees with the
    curated value, :func:`gene_disagreements` says so, and that *is* worth a curator's time.

`one_substitution`
    One residue differs from a peptide that is in the proteome, and the cell names the residue.
    **This is a description and not a defect**, which is why the verdict is spelled as the
    measurement rather than as a conclusion - "analogue" would assume engineering, and at least five
    unrelated things produce this shape, all of them correct data:

    * a **tumour neoantigen**, which is one substitution from native by definition;
    * a **germline or allelic variant**, which is the whole subject in autoimmunity;
    * a **cross-species homolog**, where a peptide conserved between human and mouse is curated under
      one species and found under the other's gene;
    * a **post-translationally modified or hybrid peptide**, where the reference proteome has the
      unmodified form;
    * an **engineered analogue**, anchor-optimised for affinity.

    There is no way to tell these apart from the sequence, and no threshold that would. The reference
    is what settles it, which is why every row carries its ``reference.id`` list: `corpus/pubmed.tsv`
    has the title and journal for 610 of them and `corpus/text_terms.tsv` the per-reference term
    counts, so the paper is one join away. A row here is an invitation to read it, not a verdict on
    it.

    **Count epitopes here, never records.** ``SLLMWITQV`` alone is 29,729 of the 36,496 records -
    **81.5 %** - so a record total is measuring one reagent choice in one antigen and not a property
    of the database. The median is 2 records per epitope and 121 of the 211 sit at one or two, which
    is what the population actually looks like.

    That largest case is also not a mislabel. ``SLLMWITQV`` is the C9V form of NY-ESO-1's
    ``SLLMWITQC``, and an anchor-optimised peptide is a **screening reagent** for the real self
    neoantigen: the ``NY-ESO-1`` label is right and those TCRs genuinely are specific for that
    antigen. What the catalogue cannot express is the *link* - the two peptides are unrelated rows -
    so a user who wants to know whether a response was found with the native peptide or a modified one
    has no way to ask. A missing cross-reference, not a defect.

    **A reference proteome is one genome and a patient cohort is not**, which is why
    :func:`by_gene` exists and is the view worth reading. One antigen yields as many distinct peptides
    as the cohort carries mutations in it, every one correctly labelled with that gene: ``PMEL`` 21
    peptides, ``INS`` 13, ``KRAS`` 12 over 50 records, ``CD20`` 8, ``p53`` 6, ``NRAS`` 3. A
    somatic-mutation panel and an antigen screen both look exactly like that, and neither is anything
    to fix. For anything tumour- or patient-derived, differing from the reference is the normal state.

`not_found`
    Neither, within one substitution. Also a question rather than a verdict: a splice junction, a
    fusion, a longer modification and a transcription error all look like this.

**Nothing here is ever rewritten, and the counts are not a defect tally.** A proteome is evidence
about a string; the publication is the authority for what the experiment used. What the table buys is
two things and no more: the three cases stop being indistinguishable, and an epitope can be named
alongside the proteome peptide it differs from, which the catalogue cannot express at all today.

Only the two self proteomes are covered. A viral or bacterial epitope's source is the pathogen's
proteome, and `antigen.species` names the pathogen there rather than the host, so strain variation is
a separate question against a per-pathogen fetch.

The output is a **committed, reviewed input** and never a build artifact, for the reason
:mod:`vdjdb.curate.promiscuity` gives: `mhcmatch` fetches the proteome from HuggingFace on first use,
so doing this inside a build would put a network call in the critical path and make the answer depend
on the day (hard rule 9). Refreshed by its own pull request through ``vdjdb antigens``.

⚠ Not to be confused with ``out/reports/epitope-sources.tsv`` (#633), which asks whether the
*corpus* is self-consistent about a peptide's source - two species or two gene labels on one string.
That one reads `chunks/` and needs no authority; this one reads a proteome and needs no corpus.
"""
from __future__ import annotations

from typing import NamedTuple

import polars as pl

#: VDJdb species -> the reference proteome `mhcmatch` fetches for it. Only the two self proteomes:
#: a viral epitope's source is the pathogen's proteome, which is a per-pathogen fetch and a different
#: question, and `antigen.species` names the *pathogen* there rather than the host.
PROTEOMES: dict[str, str] = {"HomoSapiens": "human", "MusMusculus": "mouse"}

#: How far from a proteome peptide the search goes. One, because at one substitution the hit is a
#: single identifiable peptide and the difference can be named residue by residue; at two it usually
#: is not, and a report that cannot name what differs says nothing a curator can act on.
MAX_SUBS = 1

COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "antigen.gene", "verdict",
    "source.protein", "source.gene", "source.position", "source.peptide", "source.subs",
    "records", "references",
)

#: The three verdicts. Spelled as what was measured rather than as what it might mean: a peptide one
#: residue from a proteome peptide is `one_substitution`, not "analogue", because the five things that
#: produce that shape are told apart by the reference and never by the sequence.
EXACT, ONE_SUB, NOT_FOUND = "exact", "one_substitution", "not_found"


class Source(NamedTuple):
    """One epitope against one proteome."""

    epitope: str
    species: str
    gene: str
    verdict: str
    protein: str
    source_gene: str
    position: int
    peptide: str
    subs: str
    records: int
    #: Every `reference.id` reporting this epitope, comma-joined. The row's whole purpose when the
    #: verdict is `one_substitution`: the sequence cannot say which of the five causes it is and the
    #: paper can, so the paper has to be reachable from the row.
    references: str


def sources(epitopes: pl.DataFrame, records: pl.DataFrame | None = None) -> list[Source]:
    """One row per (epitope, species) whose species has a reference proteome.

    Network-bound on first call: `mhcmatch` fetches the proteome from HuggingFace. Two batched calls
    per species - one exact, one for the leftovers at :data:`MAX_SUBS` - never one call per peptide
    (CLAUDE.md hard rule 3).
    """
    from mhcmatch.proteome import Proteome, gene_symbols
    from mhcmatch.store import fetch_proteome

    references = _references(records)
    out: list[Source] = []
    for species, organism in PROTEOMES.items():
        mine = epitopes.filter(pl.col("antigen.species") == species)
        if mine.is_empty():
            continue
        proteome = Proteome.from_hf(organism)
        # `key="name"` keys on the *whole* FASTA header, `sp|P78358|CTG1B_HUMAN`, which is exactly
        # what `SourceHit.protein` carries - so this is one lookup and no header parsing.
        symbols = gene_symbols(fetch_proteome(organism), key="name")
        peptides = mine["antigen.epitope"].unique().to_list()
        exact = proteome.find_exact_sources(peptides)
        missing = [p for p in peptides if not exact.get(p)]
        near = proteome.find_sources(missing, max_subs=MAX_SUBS, exclude_exact=True, best_only=True)

        for row in mine.iter_rows(named=True):
            peptide = row["antigen.epitope"]
            # `find_exact_sources` returns a list per peptide and `find_sources(best_only=True)` a
            # single hit or None, so both shapes arrive here and both can be empty.
            found_hit = exact.get(peptide) or near.get(peptide)
            if isinstance(found_hit, list):
                found_hit = found_hit[0] if found_hit else None
            if found_hit is None:
                out.append(Source(peptide, species, row["antigen.gene"], NOT_FOUND,
                                  "", "", -1, "", "", row["records"],
                                  references.get(peptide, "")))
                continue
            hit = found_hit
            out.append(Source(
                peptide, species, row["antigen.gene"],
                EXACT if hit.n_subs == 0 else ONE_SUB,
                hit.protein, symbols.get(hit.protein, ""), hit.position, hit.ref_peptide,
                ",".join(f"{pos + 1}{was}>{now}" for pos, now, was in hit.mutations),
                row["records"], references.get(peptide, "")))
    return out


def _references(records: pl.DataFrame | None) -> dict[str, str]:
    """``epitope -> every reference.id reporting it``, comma-joined and sorted.

    Empty when no records frame is given, which keeps the proteome lookup usable on its own; the
    column is then blank rather than absent, so the schema does not depend on the caller.
    """
    if records is None or "reference.id" not in records.columns:
        return {}
    grouped = (records.group_by("antigen.epitope")
                      .agg(pl.col("reference.id").unique().sort().str.join(",").alias("references")))
    return dict(zip(grouped["antigen.epitope"], grouped["references"], strict=True))


def table(found: list[Source]) -> pl.DataFrame:
    """The committed table: sorted so a refresh does not churn (hard rule 7)."""
    schema = {c: (pl.Int64 if c in ("source.position", "records") else pl.String) for c in COLUMNS}
    return (pl.DataFrame([tuple(s) for s in found], schema=schema, orient="row")
            .sort("antigen.species", "antigen.epitope"))


def summarise(found: pl.DataFrame) -> pl.DataFrame:
    """Per species and verdict: epitopes and the records behind them."""
    return (found.group_by("antigen.species", "verdict")
                 .agg(pl.len().alias("epitopes"), pl.col("records").sum().alias("records"))
                 .sort("antigen.species", "verdict"))


def by_gene(found: pl.DataFrame) -> pl.DataFrame:
    """How many distinct `one_substitution` peptides each curated gene carries.

    The view to read instead of a record total. A cohort carries many mutations in one antigen, so a
    gene with a dozen peptides one residue from reference is a mutation panel or an antigen screen -
    `KRAS` has 12 over 50 records - and the record count would have said almost nothing, because one
    epitope is 81.5 % of it.
    """
    return (found.filter(pl.col("verdict") == ONE_SUB)
                 .group_by("antigen.gene")
                 .agg(pl.len().alias("peptides"), pl.col("records").sum().alias("records"))
                 .sort("peptides", "records", descending=True))


def gene_disagreements(found: pl.DataFrame) -> pl.DataFrame:
    """Exactly-present epitopes whose curated gene is not the proteome's symbol.

    Case-folded, because `G6pc2` against `G6PC2` is a spelling and `lookalikes.tsv` already reports
    that class. What is left is a protein name where a gene symbol belongs, or a legacy alias - `PGT`
    for `SLCO2A1`, `IGRP` for `G6PC2` - which is what #632 asks to be able to see.
    """
    return (found.filter((pl.col("verdict") == EXACT) & (pl.col("source.gene") != "")
                         & (pl.col("antigen.gene").str.to_lowercase()
                            != pl.col("source.gene").str.to_lowercase()))
                 .select("antigen.epitope", "antigen.species", "antigen.gene", "source.gene",
                         "source.protein", "records")
                 .sort("records", descending=True))
