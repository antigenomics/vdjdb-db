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

pytestmark = pytest.mark.release

#: The epitope the validation runs on, and the motif its specific receptors are known for.
EPITOPE = "GILGFVFTL"
MOTIF = "RS"

#: A 3-mer needs this many occurrences among the epitope's documents before its lift is read. Below
#: it, one paper's handful of receptors moves the ratio and the number is noise.
MIN_OCCURRENCES = 50

#: Measured 2026-09-27 on the 2026-09 build, `documents` re-measured 2026-09-28 when `PMID_18025130`
#: landed (#161). Frozen so a change in the weighting shows up here rather than in a conclusion
#: someone draws later.
EXPECTED = {
    "documents": 662,
    "kmers_scored": 2342,
    "top_kmer": "k:IRS",
    "top_lift": 2.663,
    "median_lift": 1.170,
    "motif_kmers": 29,
    "motif_above_median": 25,
}

#: ``k:CAS`` is the germline-encoded start of nearly every beta CDR3: present in 615 of 662 documents.
#: It must not read as antigen-specific, and a weighting that made it look so would be wrong in the
#: way that matters most, because it is the first thing anyone will query.
GERMLINE_KMER = "k:CAS"
GERMLINE_TOLERANCE = 0.10


@pytest.fixture(scope="module")
def tables() -> dict[str, pl.DataFrame]:
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    needed = {n: directory / f"{n}.parquet" for n in ("records", "chains", "restriction")}
    if missing := [str(p) for p in needed.values() if not p.exists()]:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    return {n: pl.read_parquet(p) for n, p in needed.items()}


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


def test_a_germline_kmer_does_not_read_as_antigen_specific(corpus) -> None:
    """`k:CAS` starts nearly every beta CDR3, so any antigen it appears specific to is an artefact."""
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


def test_the_mhc_dictionary_recognises_all_but_the_two_known_bad_calls(tables) -> None:
    """The dictionary doubles as a proofreading report, so its findings are pinned."""
    found = mhc.unrecognised(mhc.dictionary(tables["restriction"]))
    assert sorted(found["mhc"].unique().to_list()) == ["HLA-A*08:01", "HLA-B*12"]


def test_the_corpus_is_reproducible_across_builds(tables) -> None:
    """Hard rule 7. Ids come from sorted order, so two builds of one input agree byte for byte."""
    args = (tables["records"], tables["chains"], tables["restriction"])
    first, second = build.build(*args), build.build(*args)
    for name in build.TABLES:
        assert first[name].equals(second[name]), name
