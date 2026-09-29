"""What a chunk does to the database, for the pull request that adds it.

`vdjdb qc` says whether a chunk is well formed. That is a different question from whether it is the
chunk the curator meant to submit, and the second one is answered by comparing it against the corpus
it is joining: a mistyped epitope is perfectly well formed and silently becomes a new epitope entry,
which is how `epitopes` came to hold 13 peptides under two species (#633).

So the report is written against the whole corpus rather than the file alone. Every number here needs
the rest of the database to mean anything - the score is a window over the signature columns and rises
when an independent study reports the same clonotype, and "new" means new to VDJdb, not new to the
file. That is why this reads the assembled records and not `chunks/<file>` directly.

**Nothing here fails anything.** `vdjdb qc` decides whether a chunk may land; this decides nothing and
states what the chunk does, because every finding in it has a legitimate cause as well as a careless
one. A value one character from an existing one may be a typo or may be a second stain, a second
serotype, or the same gene under two species' symbol conventions, and the curator is the only one who
knows which. A report that guessed would be wrong often enough to be ignored, and a gate that guessed
would block correct submissions.

It stops before the annotation stages. `build_master` is 9.9 s over 192,793 records; the junction,
D-posterior and segment-inference stages that follow it are another 170 s and change nothing a curator
is deciding, so the report does not wait for them.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

#: Columns whose values are compared against the rest of the corpus. A value absent from every other
#: chunk is reported - not as an error, since a genuinely new epitope is the usual reason a chunk is
#: being submitted, but because a typo looks exactly like one and nothing else would catch it.
NOVELTY_COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "antigen.gene", "mhc.a", "mhc.b", "species",
    "reference.id",
)

#: Characters removed before asking whether two values are the same value written twice. Narrow on
#: purpose: `-`, `_`, `.` and space separate nothing in a gene symbol, so `MAGE-A3` and `MAGEA3` are one
#: gene, and `PMID: 123` and `PMID:123` are one reference. **`+` and `-` are not interchangeable
#: everywhere** - measured over the corpus, a normaliser that strips every non-alphanumeric reports 11
#: collisions in `meta.cell.subset` and every one is false, because `CD8+CD27+CD45RA+CD95-` and the same
#: string ending `CD95+` are two populations. Hence a named set rather than `str.isalnum`.
_INSIGNIFICANT = str.maketrans("", "", "-_. ")


def _fold(value: str) -> str:
    """The form two spellings of one value share. Case and separators, nothing else."""
    return value.casefold().translate(_INSIGNIFICANT)


def collisions(values: list[str]) -> dict[str, list[str]]:
    """Folded form -> the two or more distinct spellings of it, largest set first.

    Edit distance was measured and is not used. Over the 427 distinct `antigen.gene` values,
    `difflib.get_close_matches` at ratio 0.85 returns 15 pairs and the great majority are real gene
    families - `MAGE-A1`/`MAGE-A2`/`MAGE-A3`/`MAGE-A4`, `PPM1`/`PPM1F`, `NBPF1`/`NBPF14`,
    `IE2`/`IE62` - while on the 72 `antigen.species` values all four pairs it finds are distinct
    organisms (`CMV`/`MCMV`/`LCMV`, `InfluenzaA`/`InfluenzaB`, `SARS-CoV`/`SARS-CoV-2`). Folding case
    and separators finds 12 genuine `antigen.gene` splits and no false one, so it is the whole check.
    """
    groups: dict[str, list[str]] = {}
    for v in values:
        if v:
            groups.setdefault(_fold(v), []).append(v)
    return {k: sorted(vs) for k, vs in groups.items() if len(vs) > 1}


#: Columns the corpus-wide look-alike scan reads. Not every column: measured over the built corpus, a
#: fold of `meta.cell.subset` reports 11 groups and all 11 are false, because a flow phenotype's `+` and
#: `-` are the information. These are the columns where a value is a *name*, so two spellings of it are
#: one thing named twice. `antigen.epitope` is included and finds nothing, which is the answer worth
#: having about the column that matters most.
LOOKALIKE_COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "antigen.gene", "mhc.a", "mhc.b", "species",
    "reference.id", "method.identification", "method.verification",
)


#: Which species column governs each scanned column. A gene symbol belongs to the **antigen's**
#: organism, not the donor's, and reading the wrong one inverts the answer: `Nef` and `NEF` sit on one
#: `species` (`HomoSapiens`, the donor) and two `antigen.species` (HIV-1 and InfluenzaA), and it is the
#: second that says whether two spellings could be two different genes.
_GOVERNING_SPECIES: dict[str, str] = {
    "antigen.epitope": "antigen.species", "antigen.gene": "antigen.species",
    "antigen.species": "antigen.species",
}


def lookalikes(records: pl.DataFrame) -> pl.DataFrame:
    """Values that differ only in case or a separator, one row per spelling. Advisory, never a gate.

    ``same.species`` is why this reports rather than decides. A fold of `antigen.gene` finds
    `MBP`/`Mbp` and `G6PC2`/`G6pc2`, and both are **correct**: HGNC capitalises human symbols and MGI
    title-cases mouse ones, so the two spellings are two species' conventions for one gene. Within one
    species there is no such excuse, which is the column a reader sorts on. It is computed against the
    species column that governs each scanned column (:data:`_GOVERNING_SPECIES`), and is ``false``
    where the corpus has no such column to check, so it never claims more than it read.
    """
    rows = []
    for col in LOOKALIKE_COLUMNS:
        if col not in records.columns:
            continue
        counts = records.group_by(col).len().sort("len", descending=True)
        n = dict(counts.iter_rows())
        sp = _GOVERNING_SPECIES.get(col, "species")
        for folded, spellings in collisions(counts[col].to_list()).items():
            one = (records.filter(pl.col(col).is_in(spellings))[sp].n_unique() <= 1
                   if sp in records.columns else False)
            for v in sorted(spellings, key=lambda v: -n[v]):
                rows.append((col, folded, v, n[v], len(spellings), one))
    return pl.DataFrame(rows, orient="row", schema={
        "column": pl.Utf8, "folded": pl.Utf8, "spelling": pl.Utf8, "records": pl.Int64,
        "spellings": pl.Int64, "same.species": pl.Boolean,
    }).sort("same.species", "records", descending=[True, True])


#: Peptides that really do occur in two organisms' proteomes, so two `(species, gene)` rows for them
#: are two independent reports rather than a defect. Named here, not re-diagnosed each run: the report
#: marks them `conserved` so a *new* collision is the only thing a reader has to look at.
#:
#: Checked one by one against the sequences, 2026-09-29 (#633). `VEALYLVCG` is in human `INS` and mouse
#: `Ins2`; `RPIIRPATL` in influenza A and B `NP`; `LRVMMLAPF` in *E. coli* and *S.* Typhi `yeiH`;
#: `LPRWYFYYL` in HCoV-HKU1 and HCoV-OC43; `KLPDDFMGC` and `TDDNALAYY` in SARS-CoV and SARS-CoV-2.
#: The set may shrink and must not grow without a sequence behind the entry.
CONSERVED_EPITOPES: frozenset[str] = frozenset({
    "VEALYLVCG", "RPIIRPATL", "LRVMMLAPF", "LPRWYFYYL", "KLPDDFMGC", "TDDNALAYY",
})


def epitope_sources(records: pl.DataFrame) -> pl.DataFrame:
    """Peptides whose source is not single-valued. Advisory, never a gate.

    Two questions, one report, because both are "which protein is this peptide from" and a reader
    wants them side by side:

    * **more than one species for one peptide.** ``epitopes`` is keyed on
      ``(antigen.epitope, antigen.species)``, so this is two rows there by design and a uniqueness
      constraint would be the wrong instrument. ``sources`` counts them.
    * **more than one gene label for one ``(peptide, species)``.** ``build_epitopes`` takes the modal
      label and its comment said ``epitopes.conflicts`` listed the rest -- a function that was never
      written, so until now the discarded labels were reported nowhere. ``genes`` counts them and
      ``antigen.gene`` is the modal one the catalogue kept.

    Three different things produce a row and only the third is a defect:

    * **a conserved peptide** - the same sequence in two organisms' proteomes, so two studies are two
      independent reports. :data:`CONSERVED_EPITOPES` names the ones checked against a sequence, and
      ``conserved`` carries it, so a reader sorts on it being false and reads what is left.
    * **a vocabulary gap** - one organism or one gene written two ways (``AdV`` beside ``HAdV5``,
      ``HEXON`` beside ``Hexon``), or ``Synthetic`` in ``antigen.species``, which is a provenance and
      not a species. These want an alias table and which way each folds is a curator's call
      (#632, #637). :func:`lookalikes` finds the case-and-separator subset of them.
    * **a mis-curation** - a peptide attributed to the wrong proteome, or a gene field holding
      something that is not a gene. ``RGPGRAFVTI`` was ``HomoSapiens`` / ``P18-I10`` on one row,
      against ``HIV-1`` / ``GP160`` on 85, where ``P18-I10`` is the laboratory name of the HIV-1
      V3-loop peptide itself; ``patches/`` answers that one now. The largest remaining is
      ``FVVKAYLPVNESFAFTADLRSNTGGQA`` with **187 gene labels**, ``Eef2`` to ``Eef188``, one per
      record - an index written into ``antigen.gene``, the same shape as the clonotype counter #625
      found in ``mhc.a``.
    """
    cols = ("antigen.epitope", "antigen.species", "antigen.gene")
    schema = {"antigen.epitope": pl.Utf8, "antigen.species": pl.Utf8, "antigen.gene": pl.Utf8,
              "sources": pl.UInt32, "genes": pl.UInt32, "records": pl.UInt32,
              "conserved": pl.Boolean}
    if any(c not in records.columns for c in cols):
        return pl.DataFrame(schema=schema)
    per_species = (
        records.group_by("antigen.epitope", "antigen.species")
        # The modal label, matching what `build_epitopes` keeps, with ties broken by sort order so
        # the two agree run to run (hard rule 7).
        .agg(pl.col("antigen.gene").mode().sort().first().alias("antigen.gene"),
             pl.col("antigen.gene").n_unique().cast(pl.UInt32).alias("genes"),
             pl.len().cast(pl.UInt32).alias("records"))
        .with_columns(pl.len().over("antigen.epitope").cast(pl.UInt32).alias("sources"))
    )
    return (per_species.filter((pl.col("sources") > 1) | (pl.col("genes") > 1))
            .with_columns(pl.col("antigen.epitope")
                          .is_in(list(CONSERVED_EPITOPES)).alias("conserved"))
            .select(*schema)
            .sort("conserved", "records", "antigen.epitope",
                  descending=[False, True, False]))


def report(files: list[str], records: pl.DataFrame) -> str:
    """Markdown: per chunk, the records and scores it contributes and the values it introduces."""
    from .anchors import noncanonical
    from .anchors import report as anchor_report

    names = sorted({Path(f).name for f in files})
    mine = records.filter(pl.col("chunk.file").is_in(names))
    rest = records.filter(~pl.col("chunk.file").is_in(names))
    out: list[str] = []

    if mine.is_empty():
        return ("No assembled records come from the changed files. Either the chunk contributed no "
                "record, or its name is not one the build reads.\n")

    out += [f"### {mine.height:,} record(s) from {len(names)} chunk(s)", "",
            "| Chunk | Records | Score 0 | Score 1 | Score 2 | Score 3 |", "|---|---:|---:|---:|---:|---:|"]
    hist = (mine.group_by("chunk.file", "vdjdb.score").len()
            .pivot(on="vdjdb.score", index="chunk.file", values="len")
            .fill_null(0).sort("chunk.file"))
    for row in hist.iter_rows(named=True):
        cells = " | ".join(str(row.get(str(s), 0)) for s in range(4))
        n = sum(row.get(str(s), 0) for s in range(4))
        out.append(f"| `{row['chunk.file']}` | {n:,} | {cells} |")
    out.append("")

    # The file's line count is not its record count: per-chunk dedup runs first.
    paired = mine.filter((pl.col("cdr3.alpha") != "") & (pl.col("cdr3.beta") != "")).height
    out += ["| Quantity | Value |", "|---|---:|",
            f"| records assembled | {mine.height:,} |",
            f"| distinct clonotypes | {mine.select('cdr3.alpha', 'cdr3.beta').n_unique():,} |",
            f"| paired (both chains) | {paired:,} |",
            f"| score 2 or 3 | {mine.filter(pl.col('vdjdb.score') >= 2).height:,} |", ""]

    out += ["### Values new to VDJdb", "",
            "A value no other chunk carries. Expected for a new epitope; a typo is indistinguishable "
            "from one, which is why they are listed.", "",
            "| Column | New values | Which |", "|---|---:|---|"]
    typos: list[str] = []
    for col in NOVELTY_COLUMNS:
        if col not in records.columns:
            continue
        known = set(rest[col].unique().to_list())
        new = sorted(v for v in mine[col].unique().to_list() if v and v not in known)
        folded = {_fold(k): k for k in known}
        typos += [f"`{col}`: **`{v}`** against `{folded[_fold(v)]}`, already in VDJdb"
                  for v in new if _fold(v) in folded]
        shown = ", ".join(f"`{v}`" for v in new[:8]) + (f" and {len(new) - 8} more" if len(new) > 8 else "")
        out.append(f"| `{col}` | {len(new)} | {shown or '-'} |")
    out.append("")

    if typos:
        out += ["#### Alert: probably a spelling of a value VDJdb already has", "",
                "**This blocks nothing.** These differ from an existing value only in case or in a "
                "`-`, `_`, `.` or space, so they will not join it: the database carries both, and every "
                "query filtering on one misses the other. `IE1` and `IE-1` are two `antigen.gene` "
                "values today for the same CMV gene.", "",
                "Two values that look alike can both be right, and only the curator knows which case "
                "this is - one stain against another, `MBP` in human against `Mbp` in mouse where the "
                "symbol conventions genuinely differ. Read it and decide; nothing here is a verdict.",
                ""]
        out += [f"* {t}" for t in typos]
        out.append("")

    # The junction anchors, against the germline of the segment each record names. No QC rule tests
    # this - they check the residue alphabet and a minimum length - and `arda.cdr3fix` repairs most of
    # them on the way through, so without this section the chunk keeps the wrong sequence and the
    # submitter never learns.
    anchor_text = anchor_report(noncanonical(mine))
    if anchor_text:
        out += [anchor_text]

    # A row keyed identically to one another chunk already reports is not a duplicate - two chunks are
    # two independent reports (CLAUDE.md) - but it is the thing that raises the score, so say so.
    key = ["cdr3.alpha", "cdr3.beta", "antigen.epitope", "mhc.a"]
    echoed = mine.join(rest.select(key).unique(), on=key, how="semi").height
    out += [f"**{echoed:,} record(s) repeat a clonotype/pMHC another chunk already reports.** That is "
            "independent replication, not duplication - it is what raises `vdjdb.score` - but a chunk "
            "where every record is an echo may be a dataset VDJdb already has.", ""]
    return "\n".join(out)
