"""One publication curated in two chunk files (#390).

`CHUNK_DEDUP_KEY` contains `reference.id`, so a group of it spanning two chunk files is one paper
reporting one clone twice - and two rows of one paper are not two independent reports, which is the
whole basis for deduplicating within a chunk. The chunk is normally the publication; where the two
come apart, the publication is what counts.
"""
from __future__ import annotations

import polars as pl

from vdjdb.io.chunks import CURATION_ONLY, _base_first, merge_repeated_references
from vdjdb.schema import CHUNK_DEDUP_KEY

CLEAN = {c: "" for c in (*CHUNK_DEDUP_KEY, "method.identification", "method.verification",
                         "meta.structure.id", "meta.epitope.id", "chunk.id", "submitter",
                         "comment")}
CLONE = {"cdr3.beta": "CASSIRSSYEQYF", "v.beta": "TRBV10-3*01", "j.beta": "TRBJ2-7*01",
         "species": "HomoSapiens", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M", "mhc.class": "MHCI",
         "antigen.epitope": "GILGFVFTL", "antigen.gene": "M", "antigen.species": "InfluenzaA",
         "reference.id": "PMID:1"}


def _frame(*rows: dict[str, str]) -> pl.DataFrame:
    return pl.DataFrame([{**CLEAN, **CLONE, **r, "chunk.row": i + 1}
                         for i, r in enumerate(rows)])


def test_one_paper_in_two_chunks_merges_and_the_submission_is_the_base():
    """The paper's own chunk is the base and the derived one fills its blanks."""
    got, report = merge_repeated_references(_frame(
        {"chunk.file": "PDB_Database.txt", "method.identification": "tetramer-sort",
         "method.verification": "structural", "meta.structure.id": "8VCX"},
        {"chunk.file": "PMID_1.txt", "method.identification": "tetramer-sort"}))
    assert got.height == 1
    assert got["chunk.file"][0] == "PMID_1.txt", "the submitted chunk is the base"
    assert got["method.verification"][0] == "structural", "the blank is filled from the other row"
    assert got["meta.structure.id"][0] == "8VCX"
    assert report["verdict"].to_list() == ["merged"]


def test_a_column_the_two_disagree_on_stops_the_merge():
    """`PDB_Database.txt` and `PMID_34433824.txt` give one clone `structural` and `tetramer-sort`:
    a solved complex and the sort that found it, which is two observations and not a duplicate."""
    got, report = merge_repeated_references(_frame(
        {"chunk.file": "PDB_Database.txt", "method.identification": "structural"},
        {"chunk.file": "PMID_1.txt", "method.identification": "tetramer-sort"}))
    assert got.height == 2, "both rows stay"
    assert report["verdict"].to_list() == ["kept"]
    assert "structural | tetramer-sort" in report["detail"][0]


def test_a_curation_column_is_not_a_disagreement():
    """`chunk.id` is a per-chunk row serial, so two chunks curating one clone always differ on it.

    Before it was excluded, all 19 corpus groups read as conflicting for that reason alone.
    """
    assert "chunk.id" in CURATION_ONLY
    got, report = merge_repeated_references(_frame(
        {"chunk.file": "PDB_Database.txt", "chunk.id": "246"},
        {"chunk.file": "PMID_1.txt", "chunk.id": "311"}))
    assert got.height == 1 and report["verdict"].to_list() == ["merged"]


def test_the_base_is_deterministic_where_neither_chunk_is_derived():
    """Two `goncharov-*` chunks are both submissions, so the tie is broken by name and never by row
    order - hard rule 7."""
    assert _base_first(["b.txt", "a.txt"]) == ["a.txt", "b.txt"]
    assert _base_first(["a.txt", "PDB_Database.txt"]) == ["a.txt", "PDB_Database.txt"]
    assert _base_first(["PDB_Database.txt", "a.txt"]) == ["a.txt", "PDB_Database.txt"]
    assert _base_first(["small_datasets_2026-05-29.txt", "zz.txt"]) == [
        "zz.txt", "small_datasets_2026-05-29.txt"]


def test_two_chunks_naming_two_papers_are_left_alone():
    """The whole point: a receptor reported by two *publications* is independent replication."""
    got, report = merge_repeated_references(_frame(
        {"chunk.file": "PMID_1.txt", "reference.id": "PMID:1"},
        {"chunk.file": "PMID_2.txt", "reference.id": "PMID:2"}))
    assert got.height == 2 and report.is_empty()


def test_the_same_clone_twice_in_one_chunk_is_not_this_pass():
    """Within-chunk duplication is `dedup`'s, and it runs first."""
    got, report = merge_repeated_references(_frame(
        {"chunk.file": "PMID_1.txt"}, {"chunk.file": "PMID_1.txt"}))
    assert got.height == 2 and report.is_empty()


def test_nothing_a_removed_row_carried_is_lost():
    """The invariant that makes the merge safe, and the one asserted on the real corpus: every
    non-blank value in the group survives on the row that is kept."""
    rows = ({"chunk.file": "PDB_Database.txt", "method.verification": "structural",
             "meta.structure.id": "6AVF"},
            {"chunk.file": "PMID_1.txt", "meta.epitope.id": "NY-ESO-160-72"})
    before = _frame(*rows)
    got, _ = merge_repeated_references(before)
    for column in before.columns:
        if column in CURATION_ONLY:
            continue
        supplied = {v for v in before[column].to_list() if v not in (None, "")}
        assert supplied <= {v for v in got[column].to_list() if v not in (None, "")}, column
