"""Every legacy proofreading failure mode, on a record built to fail it, run through both builds.

The retired build is the only specification of what VDJdb proofreading *was*, and `chunks/` cannot
demonstrate it: 231 of 231 chunks pass, so the corpus shows what both builds accept and nothing about
what either rejects. So the failures are constructed here, one record per mode, and each one is run
through :mod:`legacy_qc` - a transcription of the retired code, defects included - and through the
rules that replaced it.

Each case declares **both** verdicts and a reason for any difference. That is what makes the claim
checkable rather than asserted: a case whose new verdict changes fails here, and so does a case whose
*legacy* verdict changes, which is what stops the oracle being quietly adjusted until it agrees.

Five verdicts, and the difference between them is the whole result:

``parity``
    Both builds report the record. The chunk tier reports it under the *same name*, because the new
    rules deliberately kept the legacy finding strings, and a separate test asserts that. The master
    tier does not: legacy had one finding, ``cdr3 not C..[WF]``, where the germline-based check names
    which end is wrong and how, so the two agree that the record is broken and disagree about nothing
    else.
``stricter``
    Only the new build reports it. Either a rule with no legacy counterpart, or a row legacy waved
    through: ``species = homosapiens`` passed the retired build, which lower-cased before comparing.
``fixable``
    Legacy reported it and the new build repairs the value, so there is nothing left to report. The
    case names what the build produces and the test performs the repair; a case that claims one and
    does not perform it fails.
``legacy false positive``
    Legacy reported it and the new build is right to be silent. Only ever justified by an authority
    the case names and the test reads back - a germline, or IMGT listing an allele - never by an
    opinion.
``both silent``
    Neither build reports it, and that is deliberate on both sides rather than a gap. Murine MHC is
    the case that matters: every spelling passes on purpose.

Two coverage gates keep this exhaustive as either side grows: every legacy validator and every rule in
:data:`vdjdb.qc.rules.RULES` must be the declared finding of at least one case.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pytest

import legacy_qc as L
from vdjdb.curate import anchors
from vdjdb.curate.nomenclature import harmonise_segments, normalise_call, unresolved
from vdjdb.io.chunks import READABLE, read_chunk
from vdjdb.qc.lint import lint_file
from vdjdb.qc.rules import RULES, check
from vdjdb.schema import ALL_COLUMNS

PARITY, STRICTER, FIXABLE = "parity", "stricter", "fixable"
FALSE_POSITIVE, BOTH_SILENT = "legacy false positive", "both silent"


@dataclass(frozen=True)
class Case:
    """One deliberately broken record and what each build must say about it."""

    id: str
    row: dict[str, str]
    legacy: frozenset[str]
    new: frozenset[str]
    verdict: str
    note: str
    #: ``fixable`` only: what the new build rewrites the offending value to.
    fixed_to: str | None = None
    #: ``legacy false positive`` only: the germline residues that justify the record.
    germline: str | None = None


def _case(id: str, row: dict[str, str], legacy=(), new=(), *, verdict=PARITY, note="",
          fixed_to=None, germline=None) -> Case:
    return Case(id, row, frozenset(legacy), frozenset(new), verdict, note, fixed_to, germline)


# --------------------------------------------------------------------------------------------
# A clean record, under both builds. Every case below is this row with one field spoiled.
# --------------------------------------------------------------------------------------------

CLEAN: dict[str, str] = {
    "cdr3.beta": "CASSIRSSYEQYF", "v.beta": "TRBV10-3*01", "j.beta": "TRBJ2-7*01",
    "species": "HomoSapiens", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M", "mhc.class": "MHCI",
    "antigen.epitope": "GILGFVFTL", "antigen.gene": "M", "antigen.species": "InfluenzaA",
    "reference.id": "PMID:1",
}

# --------------------------------------------------------------------------------------------
# Tier A: py_src/ChunkQC.py, row by row
# --------------------------------------------------------------------------------------------

CHUNK_CASES: tuple[Case, ...] = (
    # --- the amino-acid alphabet and the minimum length: is_aa_seq_valid ---
    _case("cdr3-beta-has-a-residue-outside-the-twenty", {"cdr3.beta": "CASSXRSSYEQYF"},
          ["bad cdr3.beta"], ["bad cdr3.beta"],
          note="`X` for an ambiguous read. Records with unconventional residues are quarantined "
               "outside `chunks/`, so this is what arriving in the wrong directory looks like."),
    _case("cdr3-beta-is-three-residues-or-fewer", {"cdr3.beta": "CAS"},
          ["bad cdr3.beta"], ["bad cdr3.beta"],
          note="`len > 3`, so four passes and three fails. The boundary, not an arbitrary short one."),
    _case("cdr3-alpha-has-a-residue-outside-the-twenty",
          {"cdr3.alpha": "CAVRXGSQGNLIF", "v.alpha": "TRAV21*01", "j.alpha": "TRAJ42*01"},
          ["bad cdr3.alpha"], ["bad cdr3.alpha"],
          note="Both chains carry the same validator; both are covered because a rule that reads "
               "one column cannot be assumed to read the other."),
    _case("epitope-has-a-residue-outside-the-twenty", {"antigen.epitope": "GILGFVFTLX"},
          ["bad antigen.epitope"], ["bad antigen.epitope"],
          note="The epitope is checked by the same validator as a CDR3, which is why a nine-residue "
               "peptide passes the `len > 3` test that exists for junctions."),

    # --- the locus prefix on each of the five segment columns ---
    _case("v-beta-names-an-alpha-gene", {"v.beta": "TRAV21*01"},
          ["bad v.beta"], ["bad v.beta"],
          note="Chain-swapped calls are the common submission mistake, and the prefix is the only "
               "thing that catches a call put in the wrong column."),
    _case("j-beta-names-a-v-gene", {"j.beta": "TRBV10-3*01"},
          ["bad j.beta"], ["bad j.beta"], note="V pasted into the J column."),
    _case("v-alpha-names-a-beta-gene",
          {"cdr3.alpha": "CAVRDSGQGNLIF", "v.alpha": "TRBV10-3*01", "j.alpha": "TRAJ42*01"},
          ["bad v.alpha"], ["bad v.alpha"], note=""),
    _case("j-alpha-names-a-beta-gene",
          {"cdr3.alpha": "CAVRDSGQGNLIF", "v.alpha": "TRAV21*01", "j.alpha": "TRBJ2-7*01"},
          ["bad j.alpha"], ["bad j.alpha"], note=""),
    _case("d-beta-names-a-v-gene", {"d.beta": "TRBV10-3*01"},
          ["bad d.beta"], ["bad d.beta"],
          note="`d.beta` is the one segment column with no functionality rule beside it: IMGT lists "
               "two human TRBD genes and calls both functional, so the prefix is the whole check."),

    # --- the vocabularies ---
    _case("species-is-not-in-the-vocabulary", {"species": "Human"},
          ["bad species"], ["bad species"],
          note="A common name where the vocabulary wants a binomial."),
    _case("species-is-lower-case", {"species": "homosapiens"},
          [], ["bad species"], verdict=STRICTER,
          note="The retired validator was `x.lower() in speciesList`, so every casing passed and the "
               "CamelCase vocabulary was documentation rather than a rule. Downstream it is not: "
               "`IMGT_SPECIES`, the germline anchors and the recombination models are all keyed on "
               "the exact spelling, so a lower-cased row silently loses every per-species lookup."),
    _case("mhc-class-is-a-digit-one-not-a-roman-one", {"mhc.class": "MHC1"},
          ["bad mhc.class"], ["bad mhc.class"], note=""),
    _case("mhc-a-is-a-serotype-not-an-allele", {"mhc.a": "HLA-A2"},
          ["bad mhc.a"], ["bad mhc.a"],
          note="The regex constrains any value starting `HLA`, so a serological name is caught while "
               "`A2` on its own would not be."),
    _case("mhc-b-is-a-serotype-not-an-allele", {"mhc.b": "HLA-DRB1"},
          ["bad mhc.b"], ["bad mhc.b"],
          note="`mhc.b` carries its own rule. It is also where the murine class-II fragmentation "
               "lives, which is the next case."),
    _case("murine-mhc-is-unchecked-by-both", {"species": "MusMusculus", "mhc.a": "I-Ab",
                                              "mhc.b": "H2-Ab1", "mhc.class": "MHCII",
                                              "v.beta": "TRBV13-2*01", "j.beta": "TRBJ2-7*01"},
          [], [], verdict=BOTH_SILENT,
          note="Deliberate in both builds: `is_MHC_valid` returns True for anything whose first three "
               "characters are not `HLA`, which is why `I-Ab` (3,274 records), `H2-IAb` (113) and "
               "`H2-Ab1` (9) coexist. QC reads raw chunks and the corrections are declared in "
               "`patches/mhc.dict`, so a rule here would fail on the values the patch exists to fix. "
               "`assemble.epitopes.assert_mhc_resolves` gates the harmonised value instead."),
    _case("antigen-gene-is-blank", {"antigen.gene": ""},
          ["bad antigen.gene"], ["bad antigen.gene"],
          note="The only validator that is a presence check rather than a form check."),
    _case("reference-id-is-prose", {"reference.id": "see figure 3"},
          ["bad reference.id"], ["bad reference.id"],
          note="A citation has to be resolvable; `vdjdb refs` reads this column to fetch a "
               "publication year and the dashboard fails offline without one."),
    _case("reference-id-prefix-is-lower-case", {"reference.id": "pmid:1"},
          ["bad reference.id"], ["bad reference.id"],
          note="Legacy matched `PMID:` case-sensitively and `unpublished` case-insensitively. This "
               "rule carried a blanket `(?i)` and so accepted `pmid:1`, a row legacy rejected. "
               "Restored: zero corpus rows depend on the leniency."),
    _case("reference-id-says-unpublished-in-any-casing", {"reference.id": "Unpublished data"},
          [], [], verdict=BOTH_SILENT,
          note="The other half of the same rule, and the reason it cannot simply be made "
               "case-sensitive throughout."),

    # --- the three emptiness checks ---
    _case("neither-chain-has-a-cdr3", {"cdr3.beta": "", "v.beta": "", "j.beta": ""},
          ["no.cdr3"], ["no.cdr3"],
          note="The segment calls are cleared too, or `segment call with no cdr3` fires as well and "
               "the case would be testing two rules at once."),
    _case("the-epitope-is-missing", {"antigen.epitope": ""},
          ["no.antigen.seq"], ["no.antigen.seq"],
          note="A record with no epitope is not a specificity record."),
    _case("only-one-mhc-chain-is-named", {"mhc.b": ""},
          ["no.mhc"], ["no.mhc"],
          note="Both chains are required even for class I, where `mhc.b` is always B2M."),

    # --- the new rules, with no legacy counterpart ---
    _case("the-beta-junction-carries-a-second-cysteine",
          {"cdr3.beta": "CASSYCCGTEAFF"},
          [], ["internal cysteine in cdr3.beta"], verdict=STRICTER,
          note="A TCR junction has one cysteine, the Cys104 it opens with. Legacy checked the "
               "alphabet and a minimum length, so a second Cys passed. 2,804 corpus rows in 94 "
               "chunks carry one, which is why it reports and does not fail: the Jurkat receptor has "
               "one, so rare is not impossible and only the source settles it."),
    _case("the-alpha-junction-carries-a-second-cysteine",
          {"cdr3.alpha": "CAGCPRYNTDKLIF", "v.alpha": "TRAV27*01", "j.alpha": "TRAJ34*01"},
          [], ["internal cysteine in cdr3.alpha"], verdict=STRICTER,
          note="The same rule on the alpha chain, 1,521 corpus rows in 69 chunks."),

    _case("both-chains-carry-the-same-cdr3",
          {"cdr3.alpha": "CASSIRSSYEQYF", "v.alpha": "TRAV21*01", "j.alpha": "TRAJ42*01"},
          [], ["alpha and beta cdr3 identical"], verdict=STRICTER,
          note="#561. The beta sequence copied into the alpha field with the alpha calls left "
               "correct. Which chain is wrong cannot be read off the row, so it reports and does not "
               "repair."),
    _case("a-chain-has-calls-and-no-cdr3", {"v.alpha": "TRAV21*01", "j.alpha": "TRAJ42*01"},
          [], ["segment call with no cdr3"], verdict=STRICTER,
          note="The calls are information and the row is kept, but every shipped table is keyed on "
               "the CDR3, so the chain reaches no output. A submitter hears about it while the "
               "sequence can still be supplied."),
    _case("structure-id-is-a-figure-reference", {"meta.structure.id": "Fig.2"},
          [], ["structure id is not a PDB id"], verdict=STRICTER,
          note="#402. `score.confidence` awards 3 outright for a non-empty value, above every "
               "sequencing and specificity term, so a figure reference buys the top score for "
               "evidence that does not exist."),
    _case("frequency-contradicts-its-own-count",
          {"method.frequency": "0.5", "method.frequency.count": "1",
           "method.frequency.total": "10"},
          [], ["frequency disagrees with its count and total"], verdict=STRICTER,
          note="#696. `method.frequency.count` and `.total` are chunk columns now, so a submitter "
               "can report all three - and 0.5 is five times 1/10. Before #696 a cell of "
               "`method.frequency` was one shape or the other and no record could contradict "
               "itself, which is why the retired build had nothing to say here. Reports and does "
               "not repair: which of the three the paper supports is a curation question."),
    _case("method-identification-names-an-unsettled-token",
          {"method.identification": "tetramer-sort,magnetic beads"},
          [], ["undeclared method.identification token"], verdict=STRICTER,
          note="#637. The cell is a comma-separated set, so the finding is per token: "
               "`tetramer-sort` is declared and `magnetic beads` is `pending` in "
               "`proofreading/method_vocabulary.tsv`, which is what makes this row fail while a "
               "cell of only declared tokens does not. The retired build read this column as free "
               "text and said nothing about any value in it."),
    _case("v-beta-is-a-pseudogene", {"v.beta": "TRBV1*01"},
          [], ["non-functional v.beta"], verdict=STRICTER,
          note="#634. IMGT calls human TRBV1 P. Advisory rather than fatal: a P gene can rearrange, "
               "and `TRBV21-1` turns up in real repertoires on 303 chains."),
    _case("j-beta-is-an-orf", {"j.beta": "TRBJ2-2P*01"},
          [], ["non-functional j.beta"], verdict=STRICTER, note="IMGT calls TRBJ2-2P ORF."),
    _case("v-alpha-is-a-pseudogene",
          {"cdr3.alpha": "CAVRDSGQGNLIF", "v.alpha": "TRAV11*01", "j.alpha": "TRAJ42*01"},
          [], ["non-functional v.alpha"], verdict=STRICTER, note="IMGT calls human TRAV11 P."),
    _case("j-alpha-is-an-orf",
          {"cdr3.alpha": "CAVRDSGQGNLIF", "v.alpha": "TRAV21*01", "j.alpha": "TRAJ58*01"},
          [], ["non-functional j.alpha"], verdict=STRICTER,
          note="IMGT calls TRAJ58 ORF. The largest of the four by chain count, 656."),
)

# --------------------------------------------------------------------------------------------
# Tier B: runBuidDatabase.py lines 126-152, on the repaired master table
# --------------------------------------------------------------------------------------------

MASTER_CASES: tuple[Case, ...] = (
    _case("v-beta-names-a-family-and-not-a-gene", {"v.beta": "TRBV6"},
          ["gene not in IMGT"], ["gene not in IMGT"],
          note="984 chain-calls, the largest single finding. IMGT lists nine TRBV6 genes and the "
               "record picks none of them, so the row is under-specified. Only a curator can choose, "
               "which is why the new report carries the candidates rather than one of them."),
    _case("j-beta-writes-a-dot-where-a-dash-belongs", {"j.beta": "TRBJ1.2"},
          ["gene not in IMGT"], [], verdict=FIXABLE, fixed_to="TRBJ1-2",
          note="21 chains. The retired build rejected the spelling and shipped the row; the "
               "respelling rules rewrite it, so there is nothing left to report. A respelling is "
               "only applied when exactly one variant is an IMGT name for that species, which is "
               "what makes silence here safe rather than convenient."),
    _case("v-alpha-omits-the-slash-in-a-shared-gene-name", {"cdr3.alpha": "CAVRDSGQGNLIF",
                                                            "v.alpha": "TRAV29DV5",
                                                            "j.alpha": "TRAJ42*01"},
          ["gene not in IMGT"], [], verdict=FIXABLE, fixed_to="TRAV29/DV5",
          note="10 chains. `TRAV29/DV5` is one gene shared with the delta locus, so the slash is "
               "part of the name rather than a separator - which is why `_expand_slash` checks IMGT "
               "membership before reading a slash as `or`."),
    _case("the-allele-number-exceeds-the-legacy-allele-count", {"v.beta": "TRBV28*02"},
          ["allele out of range"], [], verdict=FALSE_POSITIVE, germline="TRBV28*02",
          note="`alleles_match_check` compared `int(allele)` against a per-gene count from a 741-row "
               "human immunoglobulin table, so it was a range check wearing the clothes of a "
               "membership check: `*07` passed and `*08` failed whether or not IMGT listed either. "
               "IMGT does list TRBV28*02. The one finding it produced on the whole corpus."),

    _case("the-junction-is-one-residue-short-of-the-j-anchor",
          {"cdr3.beta": "CASSIRSSYEQY"},
          ["cdr3 not C..[WF]"], ["J absent anchor"],
          note="TRBJ2-7*01 templates `SYEQYF`, so a junction ending `SYEQY` is short the Phe118. "
               "Both builds catch it; only the new one can propose the residue, because only the new "
               "one read the germline."),
    _case("the-junction-legitimately-ends-in-cysteine",
          {"cdr3.alpha": "CAVYFGNVLHC", "v.alpha": "TRAV21*01", "j.alpha": "TRAJ35*01"},
          ["cdr3 not C..[WF]"], [], verdict=FALSE_POSITIVE, germline="IGFGNVLHC",
          note="269 chain-calls. `is_qq_seq_biologically_valid` hardcoded `endswith W or F`, and "
               "TRAJ35*01 templates `IGFGNVLHC`. The junction is correct and the check was wrong."),
    _case("a-mouse-junction-legitimately-ends-in-leucine",
          {"species": "MusMusculus", "cdr3.alpha": "CAASRGSNNRLTL", "v.alpha": "TRAV6-1*01",
           "j.alpha": "TRAJ7*01", "cdr3.beta": "", "v.beta": "", "j.beta": ""},
          ["cdr3 not C..[WF]"], [], verdict=FALSE_POSITIVE, germline="DYSNNRLTL",
          note="212 chain-calls. Mouse TRAJ7*01 templates `DYSNNRLTL` where the human allele of the "
               "same number ends in Phe, so the anchor is per species and per allele and not a "
               "constant. Reading it off the reference is the only way to get this right."),
    _case("the-junction-does-not-end-in-an-anchor-and-no-j-is-named",
          {"cdr3.alpha": "CAAGGSQGNLI", "v.alpha": "TRAV29/DV5*01", "j.alpha": ""},
          ["cdr3 not C..[WF]"], ["J unanchored"],
          note="118 chains. With no J call there is no germline to read, so the check used to decline "
               "and report nothing at all - the retired build's crude rule was the only thing "
               "covering them. The universal anchor is the fallback and the finding says so, because "
               "no germline means no repair can be proposed either."),
    _case("the-junction-does-not-end-in-an-anchor-and-the-j-call-is-not-imgt",
          {"cdr3.alpha": "CAASGGYQKVT", "v.alpha": "TRAV5*01", "j.alpha": "TRAJ13-2"},
          ["cdr3 not C..[WF]", "gene not in IMGT"], ["J unanchored", "gene not in IMGT"],
          note="The same fallback reached the other way: IMGT has no `TRAJ13-2`, so there is no "
               "germline, and the call itself is the second finding."),
    _case("the-first-residue-is-not-cysteine", {"cdr3.beta": "YAISERSSYEQYF"},
          ["cdr3 not C..[WF]"], ["V corrupt anchor"],
          note="Both catch it. The new one distinguishes a mis-read Cys104 from a missing one, "
               "because no TCR folds without the disulphide: the residue is substituted rather than "
               "prepended. TRBV10-3*01 templates `CAISE`, so the three residues after the bad one "
               "place the junction and the repair is `CAISERSSYEQYF`."),
    _case("the-junction-carries-framework-behind-the-anchor",
          {"cdr3.beta": "CASSIRSSYEQYFGPG"},
          ["cdr3 not C..[WF]"], ["J under-trimmed"],
          note="Both catch it, and only one says what it is: the sequence *contains* the Phe118 and "
               "carries framework past it, so this is a submission in the wrong coordinate space "
               "rather than a truncation, and the repair trims rather than appends. Legacy saw only "
               "that the last residue was Gly."),
)


# --------------------------------------------------------------------------------------------
# Runners
# --------------------------------------------------------------------------------------------

def _write(tmp: Path, rows: list[dict[str, str]], name: str = "PMID_1.txt") -> Path:
    p = tmp / name
    # `READABLE`, not `ALL_COLUMNS`: a case may spoil an optional chunk column, which #696's
    # `method.frequency.count` is, and a header without it cannot carry the value.
    lines = ["\t".join(READABLE)]
    lines += ["\t".join(r.get(c, "") for c in READABLE) for r in rows]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _both(tmp: Path, case: Case) -> tuple[set[str], set[str]]:
    """``(legacy findings, new findings)`` for one chunk-tier case."""
    row = CLEAN | case.row
    path = _write(tmp, [row])
    legacy = {rule for _, rule in L.check_chunk([row])}
    new = set(check(read_chunk(path))["rule"].to_list())
    return legacy, new


def _master(case: Case) -> pl.DataFrame:
    """A one-row master-shaped frame for a master-tier case."""
    return pl.DataFrame([{"record_id": "VDJDB1", "chunk.file": "PMID_1.txt", "chunk.row": 1,
                          **{c: "" for c in READABLE}, **(CLEAN | case.row)}])


def _both_master(case: Case) -> tuple[set[str], set[str]]:
    """``(legacy findings, new findings)`` for one master-tier case.

    The new side is the two reports that replaced the three ``*_broken.txt`` side files:
    :func:`vdjdb.curate.nomenclature.unresolved` and :func:`vdjdb.curate.anchors.noncanonical`.
    ``unresolved`` reports per call rather than per row, so its finding is collapsed to the legacy
    name to make the two comparable.

    The legacy side reads the call **as submitted** and the new side reads it **after harmonisation**.
    That asymmetry is the measurement, not a mistake in it: ``runBuidDatabase.py`` ran ``guess_id`` and
    ``fix_both`` before its checks and no respelling at all, so a call it rejected is a call nothing in
    that build would ever have repaired. Comparing both sides post-harmonisation would make every
    ``fixable`` case read as agreement and measure nothing.
    """
    row = CLEAN | case.row
    # Harmonised first, because that is the order `build_master` runs them in and it is what makes
    # `fixable` mean anything: a call the respelling rules resolve must be gone by the time the
    # report is taken, or every fixable case would read as a finding.
    df, _ = harmonise_segments(_master(case), Path("."))
    legacy = L.check_master_row(row, Path("."))
    new: set[str] = set()
    if not unresolved(df, Path(".")).is_empty():
        new.add("gene not in IMGT")
    new |= set(anchors.noncanonical(df)["defect"].to_list())
    return legacy, new


ALL_CASES = CHUNK_CASES + MASTER_CASES


def _ids(cases: tuple[Case, ...]) -> list[str]:
    return [c.id for c in cases]


# --------------------------------------------------------------------------------------------
# Both builds report exactly what the case declares
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("case", CHUNK_CASES, ids=_ids(CHUNK_CASES))
def test_a_chunk_case_produces_the_two_declared_finding_sets(tmp_path: Path, case: Case) -> None:
    legacy, new = _both(tmp_path, case)
    assert legacy == set(case.legacy), f"legacy: {case.note}"
    assert new == set(case.new), f"new: {case.note}"


@pytest.mark.parametrize("case", MASTER_CASES, ids=_ids(MASTER_CASES))
def test_a_master_case_produces_the_two_declared_finding_sets(case: Case) -> None:
    legacy, new = _both_master(case)
    assert legacy == set(case.legacy), f"legacy: {case.note}"
    assert new == set(case.new), f"new: {case.note}"


@pytest.mark.parametrize("case", ALL_CASES, ids=_ids(ALL_CASES))
def test_the_declared_verdict_follows_from_the_two_finding_sets(case: Case) -> None:
    """The verdict is not a label. It is a relation between the two sets, and it is checked."""
    if case.verdict == PARITY:
        assert case.legacy and case.new, "parity means both builds report it"
    elif case.verdict == STRICTER:
        assert not case.legacy and case.new, "stricter means the new build alone reports it"
    elif case.verdict == BOTH_SILENT:
        assert not case.legacy and not case.new
    else:
        assert case.legacy and not case.new, f"{case.verdict} means legacy alone reported it"


CHUNK_PARITY = tuple(c for c in CHUNK_CASES if c.verdict == PARITY)


@pytest.mark.parametrize("case", CHUNK_PARITY, ids=_ids(CHUNK_PARITY))
def test_the_chunk_rules_kept_the_legacy_finding_names(case: Case) -> None:
    """`bad v.beta` still means `bad v.beta`.

    The retired build printed these strings into its warnings and a curator reading an old log has to
    find the same rule in a new report. The master tier is deliberately different: one legacy finding
    became five germline-based ones, which is the improvement rather than a rename."""
    assert case.legacy == case.new


# --------------------------------------------------------------------------------------------
# A claimed repair is performed, and a claimed false positive is justified by the reference
# --------------------------------------------------------------------------------------------

FIXABLE_CASES = tuple(c for c in ALL_CASES if c.verdict == FIXABLE)
FALSE_POSITIVES = tuple(c for c in ALL_CASES if c.verdict == FALSE_POSITIVE)


@pytest.mark.parametrize("case", FIXABLE_CASES, ids=_ids(FIXABLE_CASES))
def test_a_case_that_claims_a_repair_is_repaired(case: Case) -> None:
    """Silence is only acceptable here because the build rewrites the value. Prove it does."""
    row = CLEAN | case.row
    got = {c: normalise_call(row[c], row["species"], Path("."))
           for c in ("v.alpha", "j.alpha", "v.beta", "d.beta", "j.beta") if row.get(c)}
    assert case.fixed_to in got.values(), (
        f"{case.id} is declared fixable to {case.fixed_to!r} and the respelling rules produced "
        f"{got}")


@pytest.mark.parametrize("case", FALSE_POSITIVES, ids=_ids(FALSE_POSITIVES))
def test_a_case_that_claims_a_false_positive_names_the_reference_entry(case: Case) -> None:
    """Silence is only acceptable here because an authority says the record is right.

    Two shapes, both read back out of the reference rather than asserted: a junction is justified by
    the germline residues its J templates, and an allele by IMGT listing it.
    """
    from vdjdb.curate.nomenclature import _imgt

    row = CLEAN | case.row
    if case.germline and case.germline.startswith("TR") and "*" in case.germline:
        alleles, _ = _imgt(Path("."))[row["species"]]
        assert case.germline in alleles, f"{case.id}: IMGT does not list {case.germline}"
        return
    call = row.get("j.alpha") or row.get("j.beta")
    assert anchors.templated(row["species"], "J", call) == case.germline, (
        f"{case.id}: the germline that justifies this record is not what the case claims")
    cdr3 = row.get("cdr3.alpha") or row.get("cdr3.beta")
    assert cdr3.endswith(case.germline[-1]), (
        f"{case.id}: the junction does not end in the residue its J templates")
    assert not L.is_qq_seq_biologically_valid(cdr3), f"{case.id}: legacy would not have flagged this"


# --------------------------------------------------------------------------------------------
# Coverage: the gate that keeps "every failure mode" true as either side grows
# --------------------------------------------------------------------------------------------

def test_every_legacy_chunk_validator_has_a_case() -> None:
    declared = {f for c in ALL_CASES for f in c.legacy}
    missing = sorted(f"bad {v}" for v in L.VALIDATORS if f"bad {v}" not in declared)
    assert not missing, f"legacy validators no case exercises: {missing}"


def test_every_legacy_emptiness_and_duplicate_check_has_a_case() -> None:
    declared = {f for c in ALL_CASES for f in c.legacy} | {"duplicate"}
    for finding in ("no.cdr3", "no.antigen.seq", "no.mhc", "duplicate"):
        assert finding in declared, f"no case exercises {finding}"


def test_every_legacy_master_check_has_a_case() -> None:
    declared = {f for c in ALL_CASES for f in c.legacy}
    missing = sorted(set(L.MASTER_FINDINGS) - declared)
    assert not missing, f"legacy master-table checks no case exercises: {missing}"


#: Rules this harness structurally cannot exercise, each with where it is tested instead.
#:
#: Every case here is **one row with one field spoiled**, which is what makes a legacy verdict
#: readable beside it - the retired build also judged a row at a time. A rule that reads a *group* of
#: rows has no single-row form and no legacy counterpart to be at parity with, so a case here would be
#: a row that cannot fail it.
GROUP_RULES: dict[str, str] = {
    f"counter in {column}": "tests/unit/test_qc_advisories.py"
    for column in ("antigen.gene", "antigen.species", "mhc.a", "mhc.b")
}


def test_every_new_row_rule_has_a_case() -> None:
    """The gate that makes a new rule arrive with its own broken record, or not arrive."""
    declared = {f for c in ALL_CASES for f in c.new} | {"duplicate"} | set(GROUP_RULES)
    missing = sorted(set(RULES) - declared)
    assert not missing, (
        f"rules in vdjdb.qc.rules no case exercises: {missing}. Add a record that fails each one, "
        "with the legacy verdict beside it - or, for a rule that judges a group of rows rather than "
        f"one row, a line in GROUP_RULES naming where it is tested.")


def test_every_exempt_group_rule_exists_and_is_tested_where_it_says() -> None:
    """An exemption that names a rule nobody has, or a file nobody wrote, is a hole in the gate."""
    for rule, where in GROUP_RULES.items():
        assert rule in RULES, f"{rule} is exempt from the case gate and is not a rule"
        home = Path(where)
        assert home.exists(), f"{rule} says it is tested in {where}, which does not exist"
        column = rule.rsplit(" ", 1)[-1]
        assert column in home.read_text(), f"{where} does not mention {column}"


# --------------------------------------------------------------------------------------------
# The two checks that are a property of the file rather than of a row
# --------------------------------------------------------------------------------------------

def test_a_header_only_chunk_is_refused_by_both(tmp_path: Path) -> None:
    """`check_exist` raised `ValueError('Empty file')`. Nothing here reported it: `empty` wants a
    file with no lines at all, and the reader returns a 0-row frame without complaint."""
    path = _write(tmp_path, [])
    assert L.header_error(list(ALL_COLUMNS), 0) == "Empty file"
    assert "no-data-rows" in {f.code for f in lint_file(path)}
    assert read_chunk(path).height == 0, "the reader itself is not the gate; the lint is"


def test_a_duplicated_column_name_is_refused_by_both(tmp_path: Path) -> None:
    header = [*ALL_COLUMNS, "species"]
    path = tmp_path / "PMID_1.txt"
    row = CLEAN | {"species": "HomoSapiens"}
    path.write_text("\t".join(header) + "\n"
                    + "\t".join([*(row.get(c, "") for c in ALL_COLUMNS), "HomoSapiens"]) + "\n",
                    encoding="utf-8")
    assert L.header_error(header, 1).startswith("Duplicate columns")
    with pytest.raises(ValueError, match="duplicate columns"):
        read_chunk(path)


def test_a_missing_required_column_is_refused_by_both(tmp_path: Path) -> None:
    header = [c for c in ALL_COLUMNS if c != "method.identification"]
    path = tmp_path / "PMID_1.txt"
    path.write_text("\t".join(header) + "\n"
                    + "\t".join((CLEAN | {}).get(c, "") for c in header) + "\n", encoding="utf-8")
    assert L.header_error(header, 1).startswith("The following columns are missing")
    with pytest.raises(ValueError, match="missing required columns"):
        read_chunk(path)


def test_a_within_chunk_duplicate_is_reported_by_both(tmp_path: Path) -> None:
    """The retired *driver* never saw this one: `runBuidDatabase.py` called `drop_duplicates` before
    `ChunkQC`, so the check fired on a frame that no longer held any duplicates."""
    path = _write(tmp_path, [CLEAN, CLEAN, CLEAN])
    assert {row for row, rule in L.check_chunk([CLEAN, CLEAN, CLEAN]) if rule == "duplicate"} == {2, 3}
    found = check(read_chunk(path)).filter(pl.col("rule") == "duplicate")
    assert found["chunk.row"].to_list() == [2, 3], "the first occurrence is not a duplicate"
