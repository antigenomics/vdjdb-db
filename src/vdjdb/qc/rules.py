"""Row-level chunk validation, vectorised.

Ported from ``py_src/ChunkQC.py``, which applied each validator with ``.apply`` per cell and built
two of its three emptiness masks with ``chunk_df.T.apply`` -- transpose, then row-wise Python. Here
every rule is one polars expression over the full frame.

Two behaviours are carried over deliberately:

* **``is_MHC_valid`` passes anything not starting with ``HLA``.** The regex only constrains HLA
  spellings; murine ``H2-Kb`` and friends pass without one. That is why the murine MHC-II
  fragmentation (``I-Ab`` 3,274 vs ``H2-IAb`` 113 vs ``H2-Ab1`` 9) never tripped QC. It stays that way
  on purpose: QC reads raw chunks, and the murine and serological corrections are declared in
  ``patches/mhc.dict`` rather than edited into chunks, so a rule here would fail on the very values
  the patch exists to fix. The gate is :func:`vdjdb.assemble.epitopes.assert_mhc_resolves`, which runs
  on the harmonised values instead and fails the build on a name neither authority carries.
* **A failing chunk must exit non-zero.** The Groovy build did; the Python port replaced it with
  ``warnings.warn`` and carried on. Measured, 230 of 230 chunks pass today, so restoring the
  non-zero exit needs no quarantine list.
"""
from __future__ import annotations

import polars as pl

from ..schema import CHUNK_DEDUP_KEY, SPECIES

#: The 20 proteinogenic amino acids. Measured: all 305,031 non-empty CDR3 cells in ``chunks/`` are
#: already clean, because unconventional-residue records are quarantined outside ``chunks/``.
AA = "ARNDCQEGHILKMFPSTWYV"

_AA_SEQ = rf"^[{AA}]{{4,}}$"
#: ``HLA-<gene><digit?>*NN(:NN){0,3}``. Anything not starting with ``HLA`` is accepted -- see above.
_HLA = r"^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$"
#: The legacy prefixes, case-sensitive as legacy matched them, and ``unpublished`` anywhere in
#: the value, case-insensitive as legacy matched *that*. The blanket ``(?i)`` this carried
#: accepted ``pmid:1`` where the retired build rejected it, which is a row it caught and this
#: one waved through. Measured: zero rows in the corpus depend on the leniency.
_REFERENCE = r"^(PMID:|doi:|https?://)|(?i:unpublished)"


def _blank(col: str) -> pl.Expr:
    return pl.col(col) == ""


def _seq_ok(col: str) -> pl.Expr:
    """Empty is allowed; a present sequence must be > 3 residues of the standard 20."""
    return _blank(col) | pl.col(col).str.contains(_AA_SEQ)


def _prefix_ok(col: str, prefix: str) -> pl.Expr:
    return _blank(col) | pl.col(col).str.starts_with(prefix)


def _one_cysteine(col: str) -> pl.Expr:
    """Empty is allowed; a present junction carries no Cys after the Cys104 it opens with.

    A TCR junction has exactly one cysteine and it is the first residue. A second one is not
    impossible - the Jurkat receptor has one, and a disulphide-bonded CDR3 loop is a real thing to
    report - but it is rare enough that in ordinary submission data it reads as a transcription or
    base-calling error first. Measured over `chunks/`: 1,521 `cdr3.alpha` rows in 69 chunks and 2,804
    `cdr3.beta` rows in 94 chunks, 1.41 % of chains.

    So this is advisory: the record is kept and flagged, like the two anchor flags, and the submitter
    hears about it while the source is still to hand.
    """
    return _blank(col) | ~pl.col(col).str.slice(1).str.contains("C", literal=True)


#: ``rule id -> expression that is True when the row is GOOD``. Expressed positively so a rule
#: reads as the invariant it protects rather than as the failure it catches.
#: The chunk columns the functionality rule reads. `d.beta` is absent on purpose: IMGT lists two
#: human TRBD genes and both are functional, so the rule could never fire and a rule that cannot fire
#: is noise in every report.
SEGMENT_QC_COLUMNS: tuple[str, ...] = ("v.alpha", "j.alpha", "v.beta", "j.beta")


def _functional_ok(col: str) -> pl.Expr:
    """True unless IMGT calls the named segment ORF or P, per species. See :data:`RULES`.

    Keyed on ``species`` and the call together, because the verdict is species-specific: IMGT has
    ``TRBV7-1*01`` as ORF in human and P in *Macaca fascicularis*, which is why #634's request to
    deduplicate a tie was a cross-species artefact - keyed on ``(species, allele)`` the table has zero
    duplicate rows.
    """
    from ..curate.functionality import is_functional, verdicts

    bad = [f"{sp}\t{call}" for (sp, call), (verdict, _) in verdicts().items()
           if not is_functional(verdict)]
    key = pl.concat_str(pl.col("species"), pl.col(col), separator="\t")
    gene = pl.concat_str(pl.col("species"), pl.col(col).str.split("*").list.first(),
                         separator="\t")
    listed = [f"{sp}\t{call}" for (sp, call) in verdicts()]
    # Allele first, then the gene, and a call IMGT lists at neither depth is not this rule's finding.
    return (_blank(col)
            | (~key.is_in(bad) & key.is_in(listed))
            | (~key.is_in(listed) & ~gene.is_in(bad)))


def _method_tokens_declared() -> pl.Expr:
    """True unless ``method.identification`` names a token the vocabulary has not settled (#637).

    The specification page declares this vocabulary in prose and says "Separate phrases with a
    comma", so the unit is the **token**: 56 distinct cells over the corpus but 46 distinct tokens.
    ``proofreading/method_vocabulary.tsv`` is the list, and a token is settled when its ``status``
    is ``declared`` or ``alias``.

    **Advisory, and it must stay advisory.** A submission naming a method nobody has seen is a
    method nobody has seen, not a defect - the corpus already carries `T-Scan`, `YAMTAD system` and
    `phage display`, none of which the page ever named. What this finding buys is that the next one
    is seen when it arrives rather than counted years later, which is what #637 asked for.
    """
    from ..curate.nomenclature import method_vocabulary

    settled = [tok for tok, status, _ in method_vocabulary()
               if status in ("declared", "alias")]
    col = "method.identification"
    return (_blank(col)
            | pl.col(col).str.split(",")
                .list.eval(pl.element().str.strip_chars().is_in(settled))
                .list.all())


#: Columns that describe the **antigen or its MHC**, so a per-row counter in one of them is a
#: spreadsheet artefact rather than the column's purpose. Deliberately excludes `meta.clone.id`,
#: `meta.subject.id` and `meta.study.id`: a counter there *is* the content, and the corpus is full of
#: legitimate ones (`TCR053`..`TCR075`, `B15/S919_row007`..`row260`, donors `885`..`897`).
COUNTER_COLUMNS: tuple[str, ...] = ("antigen.gene", "antigen.species", "mhc.a", "mhc.b")

#: Distinct values a run needs before it reads as a counter rather than as two neighbouring alleles.
COUNTER_MIN_VALUES = 3

#: Share of its own integer span a run must fill. A counter is dense by construction; a gene family a
#: paper happens to report several members of is not.
COUNTER_MIN_DENSITY = 0.9


def _no_counter(col: str) -> pl.Expr:
    """True unless ``col`` holds a dense integer run inside one chunk's one epitope.

    Three confirmed instances, in three different columns, each found by hand and each costing real
    records:

    * #694, `antigen.gene`: `FVVKAYLPVNESFAFTADLRSNTGGQA` carried `Eef2`, `Eef3` ... `Eef188`, one
      label per row in lockstep with the row index. `antigen.gene` is in `CHUNK_DEDUP_KEY`, so those
      187 rows were exactly the ones that never deduplicated - 65 clonotypes held apart by a counter,
      **122 records that were never real**;
    * #625, `mhc.a`: `HLA-A*24:03` .. `HLA-A*24:20` run consecutively over 18 clonotypes and two
      epitopes, where the paper types every donor `A*24:02`. Only 4 of the 18 exist in IPD-IMGT/HLA,
      so the allele-existence gate caught 4 and the other 14 are real names carrying a wrong value;
    * the B16 chunk's `meta.epitope.id`, the same autofill frozen rather than incremented.

    Dragging a cell down a spreadsheet column increments it, so this is a recurring failure mode and
    not three accidents. What makes it detectable without judgement is that these columns are
    properties of the peptide: within one paper and one epitope they should be constant, and a dense
    run of `prefix`+integer is not a curator reporting two alleles.

    Advisory. The one current finding is #625's, already declared in `patches/mhc.dict` and repaired
    at build time, so the rule starts as a regression guard on a clean corpus - which is the state
    #597 says a new rule should start from.
    """
    group = ["chunk.file", "antigen.epitope"]
    prefix = pl.col(col).str.extract(r"^(.*?)\d+$", 1)
    number = pl.col(col).str.extract(r"(\d+)$", 1).cast(pl.Int64, strict=False)
    span = number.max().over(group) - number.min().over(group) + 1
    return ~(
        (pl.col(col).n_unique().over(group) >= COUNTER_MIN_VALUES)
        & (prefix.n_unique().over(group) == 1)
        & prefix.is_not_null()
        & ((number.n_unique().over(group) / span) >= COUNTER_MIN_DENSITY)
    )


RULES: dict[str, pl.Expr] = {
    "bad cdr3.alpha": _seq_ok("cdr3.alpha"),
    "bad cdr3.beta": _seq_ok("cdr3.beta"),
    "bad antigen.epitope": _seq_ok("antigen.epitope"),
    "bad v.alpha": _prefix_ok("v.alpha", "TRAV"),
    "bad j.alpha": _prefix_ok("j.alpha", "TRAJ"),
    "bad v.beta": _prefix_ok("v.beta", "TRBV"),
    "bad d.beta": _prefix_ok("d.beta", "TRBD"),
    "bad j.beta": _prefix_ok("j.beta", "TRBJ"),
    "bad species": pl.col("species").is_in(list(SPECIES)),
    "bad mhc.a": _blank("mhc.a") | ~pl.col("mhc.a").str.starts_with("HLA")
                 | pl.col("mhc.a").str.contains(_HLA),
    "bad mhc.b": _blank("mhc.b") | ~pl.col("mhc.b").str.starts_with("HLA")
                 | pl.col("mhc.b").str.contains(_HLA),
    "bad mhc.class": pl.col("mhc.class").is_in(["MHCI", "MHCII"]),
    "bad antigen.gene": ~_blank("antigen.gene"),
    "bad reference.id": _blank("reference.id") | pl.col("reference.id").str.contains(_REFERENCE),
    # An internal cysteine. The two anchor flags `v.canonical`/`j.canonical` already mark a junction
    # that does not open with Cys104 or close with Phe/Trp118; neither notices a Cys in the middle,
    # and no shipped column did until this rule. Advisory, and reported at the import stage, which is
    # where the submitter can still check it against their own source.
    "internal cysteine in cdr3.alpha": _one_cysteine("cdr3.alpha"),
    "internal cysteine in cdr3.beta": _one_cysteine("cdr3.beta"),
    "no.cdr3": ~(_blank("cdr3.alpha") & _blank("cdr3.beta")),
    "no.antigen.seq": ~_blank("antigen.epitope"),
    "no.mhc": ~(_blank("mhc.a") | _blank("mhc.b")),
    # #561. A paired record whose two chains have the same CDR3 is a transcription error: the
    # beta sequence copied into the alpha field, with the V and J calls left correct. Which chain is
    # wrong cannot be known from the row, so this reports and does not repair -- 99 records on the
    # current corpus, 98 of them from two references. Advisory, so the build does not fail on a
    # defect only a curator can fix.
    "alpha and beta cdr3 identical": (_blank("cdr3.alpha") | _blank("cdr3.beta")
                                      | (pl.col("cdr3.alpha") != pl.col("cdr3.beta"))),
    # A V or J call for a chain whose CDR3 is blank. The call is real information and the row is
    # kept, but the chain cannot reach any output: every shipped table is keyed on the CDR3, so the
    # segment is carried in `chunks/` and dropped by the build, which is the kind of silent loss a
    # submitter should hear about while they can still fix it.
    #
    # Two ways a row gets here. A submission that names the genes and omits the sequence, which is
    # the case this flags for the submitter. And #561, where the alpha CDR3 was a copy of the beta
    # and was cleared, leaving the paper's genuine alpha V and J behind on purpose -- that is a
    # deliberate state, recorded in the row's `comment`, and the flag keeps it visible rather than
    # letting it look like a clean record.
    #
    # Advisory: it reports a chain that cannot ship, not a row that is wrong.
    "segment call with no cdr3": ~(
        (_blank("cdr3.alpha") & ~(_blank("v.alpha") & _blank("j.alpha")))
        | (_blank("cdr3.beta") & ~(_blank("v.beta") & _blank("j.beta")))
    ),

    # `meta.structure.id` is specified as "PDB structure ID if one exists, blank otherwise"
    # (docs/standards/chunk-format.md), and `score.confidence` reads it as the strongest evidence
    # there is: a non-empty value awards 3 outright, above every sequencing and specificity term,
    # because a solved TCR:pMHC complex is direct proof of binding.
    #
    # So the field is not free text, and anything in it that is not a PDB id awards the top score
    # for evidence that does not exist. Measured 2026-09-28: 2,765 of the 3,246 rows carrying a
    # value hold a figure or table reference instead - `Fig.2, Fig. 3, Fig.4, ...` 2,352 rows,
    # `Fig 9, Supp Fig 5, Supp Table 5-8` 400, `Fig3b,Fig3c` 12, `56I` 1. Because the score is a
    # maximum over the sample signature, they reach much further than their own rows: **4,817 of
    # the 6,007 records that score 3 owe it to one of those values**, and blanking the field would
    # take 3,752 of them to 0 and 1,065 to 1.
    #
    # The shape is the PDB entry id: a digit then three alphanumerics, case-insensitive. That is the
    # whole check; whether the entry exists and contains the receptor is a network question and
    # belongs in the structure pass of #402, not in a chunk gate that must run offline.
    #
    # Advisory: the corpus fails it on 2,765 rows today, and what to do with those is a curation
    # decision rather than a submission error.
    "structure id is not a PDB id": (
        _blank("meta.structure.id")
        | pl.col("meta.structure.id").str.contains(r"^[0-9][A-Za-z0-9]{3}$")
    ),
    # #634. IMGT's own F / ORF / P verdict on the segment the row names, resolved allele-first and
    # then gene-level, per species. No other rule asks this: the two above it ask whether the call
    # *looks* like a TRBV name, and `curate.nomenclature` asks whether IMGT *has* it.
    #
    # Advisory, and that is about the biology rather than caution. A pseudogene V call is not
    # automatically wrong - a P gene can rearrange, and `TRBV21-1` (303 chains) turns up in real
    # repertoires; `annotate/junction.py` already names it as a pseudogene the recombination model
    # has no allele for. IMGT also reclassifies genes between releases, so a gate would fail on a
    # reference update rather than on a curation error. Measured 2026-09-29 over the built corpus:
    # 2,608 chain-segments, 1,172 V and 1,436 J, the largest being `TRAJ58*01` ORF on 656 chains.
    # `out/reports/functionality.tsv` is the per-chain form with IMGT's spelling and whether the
    # verdict came from the allele or from its gene.
    **{f"non-functional {col}": _functional_ok(col) for col in SEGMENT_QC_COLUMNS},
    # A spreadsheet counter in a column that describes the antigen (#694, #625). See `_no_counter`
    # for the three confirmed instances and why the `meta.*` identifier columns are excluded.
    **{f"counter in {col}": _no_counter(col) for col in COUNTER_COLUMNS},
    # A `method.identification` token the vocabulary has not settled (#637). Advisory, deliberately.
    "undeclared method.identification token": _method_tokens_declared(),
}


def check(df: pl.DataFrame) -> pl.DataFrame:
    """One row per (chunk file, chunk row, failing rule).

    ``duplicate`` marks the second and later rows sharing a :data:`CHUNK_DEDUP_KEY` within one
    chunk. Those are removed by the reader's per-chunk deduplication, so they are reported as a
    curation signal rather than an error that blocks the build.
    """
    flagged = df.with_columns(
        *(expr.not_().alias(f"__{rule}") for rule, expr in RULES.items()),
        (pl.int_range(pl.len()).over(["chunk.file", *CHUNK_DEDUP_KEY]) > 0).alias("__duplicate"),
    )
    cols = [f"__{r}" for r in RULES] + ["__duplicate"]
    return (
        flagged.select("chunk.file", "chunk.row", *cols)
        .unpivot(index=["chunk.file", "chunk.row"], on=cols,
                 variable_name="rule", value_name="failed")
        .filter("failed")
        .with_columns(pl.col("rule").str.strip_prefix("__"))
        .drop("failed")
        .sort("chunk.file", "chunk.row", "rule")
    )


def summarise(findings: pl.DataFrame) -> pl.DataFrame:
    """Findings per rule, most frequent first."""
    return (findings.group_by("rule")
            .agg(pl.len().alias("rows"), pl.col("chunk.file").n_unique().alias("chunks"))
            .sort("rows", descending=True))
