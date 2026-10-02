"""The corpus against the real database: does the instrument recover something already known?

Marked ``release``: needs a built directory (``VDJDB_TABLES``, default ``out/tables``). The unit tests
check that the weighting is arithmetically right, against ``sklearn``. These check that it is
*useful*, which arithmetic cannot: a tf-idf corpus over receptor k-mers either recovers the motifs
immunology already documents for an epitope, or it is a table nobody should draw a conclusion from.

The case is GILGFVFTL, the influenza A M1 epitope, whose specific TCRs are known for an RS motif in
the beta CDR3. Nothing in the build knows that, so recovering it is evidence the lift means something.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.corpus import build, mhc, query, tokens
from vdjdb.emit.vdjdb3 import read_table

pytestmark = pytest.mark.release

#: The epitope the validation runs on, and the motif its specific receptors are known for.
EPITOPE = "GILGFVFTL"
MOTIF = "RS"

#: A 3-mer needs this many occurrences among the epitope's documents before its lift is read. Below
#: it, one paper's handful of receptors moves the ratio and the number is noise.
MIN_OCCURRENCES = 50

#: Measured 2026-09-27 on the 2026-09 build, `documents` re-measured 2026-09-28 when `PMID_18025130`
#: landed (#161) and again 2026-09-29 when `PMID: 34433824` lost its space (#637): the two spellings
#: were two documents for one paper, sharing the epitope `GQVELGGGNAVEVCK`, so the count falls by one
#: rather than rising. Frozen so a change in the weighting shows up here rather than in a conclusion
#: someone draws later.
#:
#: `motif_above_median` re-measured 2026-09-30 on `annotate_junctions`: **25 -> 26**. The corpus's
#: `v:` and `j:` token families are the shipped segment calls, and arda 2.36 re-calls a segment whose
#: germline the junction contradicts, so the document-frequency weighting moves. The direction is the
#: one to want - one more of the 29 `RS` k-mers sits above the median, so the family's claim is
#: stronger rather than weaker - and `top_kmer`, `top_lift` and `median_lift` did not move at all.
#:
#: Re-measured 2026-10-01 on #390: **26 -> 25**, a declared trade rather than a drift. Merging the
#: 18 rows where one publication was curated in two chunk files changes the document frequency of
#: the tokens those rows carried, and one of the 29 `RS` k-mers crosses back below the median. Every
#: other measurement here is unchanged - `documents` 661, `kmers_scored` 2,342, `top_kmer` `k:IRS`,
#: `top_lift` 2.663, `median_lift` 1.170 - and the claim the test exists for is the *family* sitting
#: above the middle, which at 25 of 29 is 0.862 against the 0.8 floor the last assertion pins.
#: A duplicate row is not evidence, so removing it is the right answer even where a derived
#: statistic reads marginally weaker for it.
#: PDB primary-citation reconciliation (#839): 23 structure references become nine new
#: paper references, reducing documents from 661 to 647 without removing any of 391 PDB rows.
#: The other expectations remain unchanged; the full CI run passed their assertions.
#: #845 restores metadata-distinct observations: 26 of 29 RS k-mers are above
#: the median (previously 25 of 29); the other acceptance expectations still pass.
EXPECTED = {
    "documents": 646,  # PMID:41286513 adds one publication after the PDB corrections.
    "kmers_scored": 2342,
    "top_kmer": "k:IRS",
    "top_lift": 2.663,
    "median_lift": 1.170,
    "motif_kmers": 29,
    "motif_above_median": 26,
}

#: ``k:CAS`` is the germline-encoded start of nearly every beta CDR3: present in 614 of 661 documents.
#: It must not read as enriched under any condition, and a weighting that made it look so would be
#: wrong in the way that matters most, because it is the first thing anyone will query.
GERMLINE_KMER = "k:CAS"
GERMLINE_TOLERANCE = 0.10


@pytest.fixture(scope="module")
def tables() -> dict[str, pl.DataFrame]:
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    needed = {n: directory / f"{n}.parquet" for n in ("records", "chains", "restriction")}
    if missing := [str(p) for p in needed.values() if not p.exists()]:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    return {n: read_table(directory, n) for n in needed}


@pytest.fixture(scope="module")
def corpus(tables: dict[str, pl.DataFrame]) -> dict[str, pl.DataFrame]:
    return build.build(tables["records"], tables["chains"], tables["restriction"])


@pytest.fixture(scope="module")
def scored(corpus: dict[str, pl.DataFrame]) -> list[tuple[float, str]]:
    """Every CDR3 3-mer's lift on the epitope's documents, above the occurrence floor, best first.

    One `lift_family` call rather than a loop over `lift`. The loop recomputed the condition group and
    the family totals once per term - the same values every time - and cost **50.7 s**, a quarter of
    the whole suite, against **0.01 s** batched (`ROADMAP_local.md` §57.1). The two agree term for
    term, which `tests/unit/test_corpus_query.py` asserts.
    """
    scored = query.lift_family(corpus, "cdr3_kmer", [f"e:{EPITOPE}"], over="occurrences",
                               min_units=MIN_OCCURRENCES)
    return [(r["lift"], r["term"]) for r in scored.iter_rows(named=True)
            if r["lift"] is not None]


def test_the_corpus_covers_every_offline_family(corpus) -> None:
    """The text family needs a network refresh; the other eleven come from `chunks/`."""
    present = set(corpus["terms"]["family"].to_list())
    offline = {name for name in tokens.FAMILIES if name != "word"}
    assert offline <= present, f"missing: {sorted(offline - present)}"
    assert corpus["documents"].height == EXPECTED["documents"]


def test_every_document_is_l2_normalised(corpus) -> None:
    norms = (corpus["postings"].group_by("document_id")
                               .agg((pl.col("weight") ** 2).sum().sqrt().alias("norm")))
    assert norms["norm"].min() == pytest.approx(1.0, abs=1e-9)
    assert norms["norm"].max() == pytest.approx(1.0, abs=1e-9)


def test_the_known_motif_is_the_highest_lifting_kmer_for_its_epitope(scored) -> None:
    """The acceptance criterion. Nothing in the build knows GILGFVFTL has an RS motif.

    Of the 2,342 CDR3 3-mers with 50 or more occurrences among this epitope's documents, the top one is
    `k:IRS` at 2.66x against a median of 1.17x.
    """
    assert len(scored) == EXPECTED["kmers_scored"]
    top_lift, top_term = scored[0]
    assert top_term == EXPECTED["top_kmer"]
    assert top_lift == pytest.approx(EXPECTED["top_lift"], abs=0.01)
    assert MOTIF in top_term, "the top k-mer must carry the motif, not merely be reproducible"


def test_the_motif_family_sits_above_the_middle_of_the_distribution(scored) -> None:
    """One high k-mer could be chance. The whole RS family being high is the claim."""
    lifts = sorted(v for v, _ in scored)
    median = lifts[len(lifts) // 2]
    assert median == pytest.approx(EXPECTED["median_lift"], abs=0.01)
    motif = [(v, t) for v, t in scored if MOTIF in t[len("k:"):]]
    assert len(motif) == EXPECTED["motif_kmers"]
    above = sum(1 for v, _ in motif if v > median)
    assert above == EXPECTED["motif_above_median"]
    assert above / len(motif) > 0.8, "a motif family scattered around the median is not a signal"


def test_a_germline_kmer_reads_flat_under_every_kind_of_condition(corpus) -> None:
    """`k:CAS` starts nearly every beta CDR3, so any condition it appears enriched under is an artefact.

    Both kinds of condition are asserted, and they are not the same question. `e:GILGFVFTL` names an
    epitope the receptors were actually shown; `a:HIV-1` is provenance, a union over every epitope of
    that species and every restriction, which is a legitimate axis but not a specificity one. A
    germline k-mer has to read flat under both, which is a stronger claim than either alone.
    """
    for condition in ("a:HIV-1", f"e:{EPITOPE}"):
        got = query.lift(corpus, GERMLINE_KMER, [condition], over="occurrences")
        assert got.lift is not None
        assert abs(got.lift - 1.0) < GERMLINE_TOLERANCE, str(got)


def test_a_document_level_lift_cannot_answer_a_common_kmer(corpus) -> None:
    """Why the occurrence mode exists, asserted rather than only written down.

    `k:CAS` is in 614 of 661 documents, so its document-level lift is bounded near 1 no matter how
    specific it is. A caller reaching for the document mode on a common token gets a number with no
    room to move, and this test is what says so.
    """
    by_document = query.lift(corpus, GERMLINE_KMER, [f"e:{EPITOPE}"], over="documents")
    assert by_document.term_units / by_document.units > 0.9
    assert by_document.lift is not None
    assert by_document.lift < 1.1


def test_a_search_returns_the_papers_that_report_the_epitope(corpus, tables) -> None:
    """The endpoint's own question, and what the epitope k-mer family buys on top of it.

    Nine of the top ten report `GILGFVFTL` exactly. The tenth, `PMID:27036003`, reports `GILEFVFTL`
    and `GILGLVFTL`, single-residue variants of it, and ranks because the `ek:` family matches their
    shared 3-mers. An exact-epitope search cannot find an altered-peptide-ligand study of the epitope
    it is asking about; this is the reason `ek:` is in the vocabulary rather than `e:` alone.
    """
    hits = query.score(corpus, query.refsearch_query(epitope=EPITOPE), limit=10)
    assert hits.height == 10
    records = tables["records"]
    exact = set(records.filter(pl.col("antigen.epitope") == EPITOPE)["reference.id"])
    ranked = hits["reference.id"].to_list()
    assert sum(1 for r in ranked if r in exact) >= 9

    windows = {EPITOPE[i:i + tokens.K] for i in range(len(EPITOPE) - tokens.K + 1)}
    for reference in (r for r in ranked if r not in exact):
        epitopes = records.filter(pl.col("reference.id") == reference)["antigen.epitope"].to_list()
        assert any(e[i:i + tokens.K] in windows
                   for e in epitopes for i in range(len(e) - tokens.K + 1)), (
            f"{reference} ranks for {EPITOPE} while sharing no {tokens.K}-mer with it")


def test_the_mhc_dictionary_recognises_every_call(tables) -> None:
    """The dictionary doubles as a proofreading report, so its findings are pinned - at none.

    It used to pin two: `HLA-A*08:01`, for which no `HLA-A*08` exists at any resolution, and
    `HLA-B*12`, a serological broad antigen that split into B*44 and B*45. Both are corrected in
    `patches/mhc.dict`, and `assert_mhc_resolves` now fails the build on an unresolved call, so this
    list can only be empty on a build that completed. Empty is the useful state: a row here means a
    call that arrived since, not one of two strings a reader has to remember to ignore.
    """
    found = mhc.unrecognised(mhc.dictionary(tables["restriction"]))
    assert found["mhc"].unique().to_list() == []


def test_the_corpus_is_reproducible_across_builds(tables) -> None:
    """Hard rule 7. Ids come from sorted order, so two builds of one input agree byte for byte."""
    args = (tables["records"], tables["chains"], tables["restriction"])
    first, second = build.build(*args), build.build(*args)
    for name in build.TABLES:
        assert first[name].equals(second[name]), name
