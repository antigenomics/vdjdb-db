"""The derived id levels: clonotype, clone, pMHC and epitope.

``ROADMAP.md`` section 10.5 names seven invariants. Six of them are checks over a built directory and
are exercised through :mod:`vdjdb.identity.checks` below. The seventh is a property of two builds, so
it is the one test here that does the assembly twice: reading the chunks in a permuted order, and
adding a chunk and removing it again, must leave every other id untouched.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.identity import checks, lifecycle
from vdjdb.identity.levels import (
    CLONE,
    CLONOTYPE,
    EPITOPE,
    ID_HEX,
    LEVELS,
    PMHC,
    clone_ids,
    derive,
    duplicate_ids,
    hash_id,
    join_key,
    level_of,
    recomputed_mismatches,
)

RECORD = {"record_id": "VDJDB0000000001", "species": "HomoSapiens", "mhc.a": "HLA-A*02:01",
          "mhc.b": "B2M", "antigen.epitope": "GILGFVFTL", "reference.id": "PMID:1"}
CHAIN = {"record_id": "VDJDB0000000001", "gene": "TRB", "cdr3": "CASSIRSSYEQYF",
         "v.segm": "TRBV10-3*01", "j.segm": "TRBJ2-7*01", "TCR_hash": "abc",
         # `species` is on `records` in the shipped tables; it rides on the fixture so `derive` can
         # be called directly, and the `checks` tests below drop it to exercise the join path.
         "species": "HomoSapiens"}


def records(*overrides: dict) -> pl.DataFrame:
    return pl.DataFrame([RECORD | o for o in (overrides or ({},))])


def chains(*overrides: dict) -> pl.DataFrame:
    return pl.DataFrame([CHAIN | o for o in (overrides or ({},))])


# --------------------------------------------------------------------------------------------
# The id itself
# --------------------------------------------------------------------------------------------

def test_an_id_is_sha256_of_its_key_and_not_a_dependency_hash() -> None:
    """Frozen against a literal, so a polars or Python upgrade cannot renumber the database.

    ``pl.Expr.hash`` is xxhash and polars does not specify its output across versions. This id ships
    to consumers, so the algorithm has to be one we own (``ROADMAP.md`` section 10.3).
    """
    import hashlib

    key = join_key(["HomoSapiens", "TRB", "CASSIRSSYEQYF", "TRBV10-3*01", "TRBJ2-7*01"])
    assert key == "HomoSapiens\x1fTRB\x1fCASSIRSSYEQYF\x1fTRBV10-3*01\x1fTRBJ2-7*01"
    want = "CT" + hashlib.sha256(key.encode()).hexdigest()[:ID_HEX]
    assert hash_id("CT", key) == want
    assert hash_id("CT", key) == "CT2df51ea0980263d2"


def test_every_prefix_is_distinct_and_resolves_to_its_level() -> None:
    prefixes = [lv.prefix for lv in LEVELS]
    assert len(set(prefixes)) == len(prefixes)
    assert level_of("CT3efa364e7f84aac9") == "clonotype"
    assert level_of("CX0000b5a8beefa2ae") == "clone"
    assert level_of("PM22b0715779e2a79a") == "pmhc"
    assert level_of("EPc9eb36fceaab0098") == "epitope"
    assert level_of("VDJDB0000000001") == "record"
    assert level_of("ZZdeadbeef") is None


def test_a_missing_key_column_raises_rather_than_hashing_a_column_of_blanks() -> None:
    """Otherwise one id would cover the whole table and every downstream join would still work."""
    with pytest.raises(KeyError, match="clonotype"):
        derive(chains().drop("v.segm"), CLONOTYPE)


# --------------------------------------------------------------------------------------------
# What changes an id and what does not
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("level,frame_fn,change", [
    (CLONOTYPE, chains, {"cdr3": "CASSIRSSYEQYFF"}),
    (CLONOTYPE, chains, {"v.segm": "TRBV10-3*02"}),
    (CLONOTYPE, chains, {"gene": "TRA"}),
    (PMHC, records, {"mhc.a": "HLA-A*02:02"}),
    (PMHC, records, {"antigen.epitope": "NLVPMVATV"}),
    (EPITOPE, records, {"antigen.epitope": "NLVPMVATV"}),
])
def test_a_key_change_produces_a_different_id(level, frame_fn, change) -> None:
    out = derive(frame_fn({}, change), level)
    assert out[level.column][0] != out[level.column][1]


@pytest.mark.parametrize("level,frame_fn,change", [
    (CLONOTYPE, chains, {"TCR_hash": "different"}),
    (CLONOTYPE, chains, {"record_id": "VDJDB0000000002"}),
    (PMHC, records, {"reference.id": "PMID:2"}),
    (EPITOPE, records, {"mhc.a": "HLA-B*07:02"}),
])
def test_a_non_key_change_keeps_the_id(level, frame_fn, change) -> None:
    """Re-annotating a record is not a new clonotype, and a second paper is not a second epitope."""
    out = derive(frame_fn({}, change), level)
    assert out[level.column][0] == out[level.column][1]


def test_the_pmhc_key_ignores_mhc_class_because_the_alleles_determine_it() -> None:
    out = derive(records({}, {"mhc.class": "MHCII"}), PMHC)
    assert out[PMHC.column][0] == out[PMHC.column][1]


# --------------------------------------------------------------------------------------------
# The clone
# --------------------------------------------------------------------------------------------

def _paired(record_id: str, alpha: str, beta: str) -> pl.DataFrame:
    return pl.DataFrame([
        CHAIN | {"record_id": record_id, "gene": "TRA", "clonotype_id": alpha},
        CHAIN | {"record_id": record_id, "gene": "TRB", "clonotype_id": beta},
    ])


def test_a_clone_id_does_not_depend_on_which_chain_came_first() -> None:
    """The pair is sorted before hashing, so a reader meeting beta first gets the same clone."""
    a = clone_ids(_paired("VDJDB0000000001", "CTaaa", "CTbbb"))
    b = clone_ids(_paired("VDJDB0000000001", "CTbbb", "CTaaa"))
    assert a["clone_id"][0] == b["clone_id"][0]


def test_two_papers_reporting_one_pair_report_one_clone() -> None:
    both = pl.concat([_paired("VDJDB0000000001", "CTaaa", "CTbbb"),
                      _paired("VDJDB0000000002", "CTaaa", "CTbbb")])
    out = clone_ids(both)
    assert out.height == 2
    assert out["clone_id"].n_unique() == 1


def test_a_single_chain_record_gets_no_clone_row() -> None:
    """Having no clone is the curated state; `build_chains` fills the column with an empty string."""
    one = chains({"clonotype_id": "CTaaa"})
    assert clone_ids(one).is_empty()


def test_two_chains_of_the_same_gene_are_not_a_clone() -> None:
    """A record with two betas is a curation defect, not a pair, so it gets no clone id."""
    two_betas = pl.DataFrame([
        CHAIN | {"gene": "TRB", "clonotype_id": "CTaaa"},
        CHAIN | {"gene": "TRB", "clonotype_id": "CTbbb"},
    ])
    assert clone_ids(two_betas).is_empty()


# --------------------------------------------------------------------------------------------
# Invariant 3: order independence
# --------------------------------------------------------------------------------------------

def _rows(n: int) -> list[dict]:
    return [CHAIN | {"record_id": f"VDJDB{i:010d}", "cdr3": f"CASS{'A' * (i % 5)}F",
                     "v.segm": f"TRBV{i % 3 + 1}*01"} for i in range(n)]


def test_a_permuted_read_order_changes_no_id() -> None:
    """Invariant 3, the property the chunk-order question asks for.

    A derived id is a hash of its own key, so shuffling the rows can only shuffle the output. A
    counter would renumber every row after the first difference, which is why none of these levels
    uses one.
    """
    import random

    rows = _rows(40)
    forward = derive(pl.DataFrame(rows), CLONOTYPE)
    shuffled = rows[:]
    random.Random(20260927).shuffle(shuffled)
    reverse = derive(pl.DataFrame(shuffled), CLONOTYPE)
    a = dict(zip(forward["record_id"], forward[CLONOTYPE.column], strict=True))
    b = dict(zip(reverse["record_id"], reverse[CLONOTYPE.column], strict=True))
    assert a == b


def test_adding_and_removing_a_chunk_changes_no_other_id() -> None:
    """Invariant 3, the second half. Also the reason a registry is not consulted for these levels."""
    base = _rows(20)
    extra = [CHAIN | {"record_id": "VDJDB0000000099", "cdr3": "CAWSVNEW", "v.segm": "TRBV9*01"}]
    before = derive(pl.DataFrame(base), CLONOTYPE)
    during = derive(pl.DataFrame(base + extra), CLONOTYPE)
    after = derive(pl.DataFrame(base), CLONOTYPE)

    def by_record(frame: pl.DataFrame) -> dict[str, str]:
        return dict(zip(frame["record_id"], frame[CLONOTYPE.column], strict=True))

    assert by_record(before) == by_record(after)
    assert by_record(before).items() <= by_record(during).items()
    assert "VDJDB0000000099" in by_record(during)


# --------------------------------------------------------------------------------------------
# Invariants 1 and 2 catch a planted fault
# --------------------------------------------------------------------------------------------

def test_a_collision_is_reported_rather_than_joined_through() -> None:
    """Invariant 1. Two distinct keys carrying one id is what a truncated digest would look like."""
    out = derive(chains({}, {"cdr3": "CASSDIFFERENT"}), CLONOTYPE)
    forged = out.with_columns(pl.lit("CTcollision").alias(CLONOTYPE.column))
    assert duplicate_ids(forged, CLONOTYPE.column, CLONOTYPE.key).height == 1
    assert duplicate_ids(out, CLONOTYPE.column, CLONOTYPE.key).is_empty()


def test_an_id_that_is_not_its_own_key_is_reported() -> None:
    """Invariant 2. A bad join attaches a plausible id, and only recomputing catches it."""
    out = derive(chains(), CLONOTYPE)
    assert recomputed_mismatches(out, CLONOTYPE).is_empty()
    wrong = out.with_columns(pl.lit("CT0000000000000000").alias(CLONOTYPE.column))
    assert recomputed_mismatches(wrong, CLONOTYPE).height == 1


def test_check_passes_on_a_consistent_pair_of_tables() -> None:
    rec = derive(derive(records(), PMHC), EPITOPE)
    ch = derive(chains(), CLONOTYPE).with_columns(pl.lit("").alias(CLONE.column)).drop("species")
    assert checks.check({"records": rec, "chains": ch}) == []


def test_check_reports_a_chain_whose_record_does_not_exist() -> None:
    rec = derive(derive(records(), PMHC), EPITOPE)
    ch = derive(chains({"record_id": "VDJDB0000000404"}), CLONOTYPE).drop("species")
    found = checks.check({"records": rec, "chains": ch})
    assert any(f.invariant == 6 and "record_id that records does not have" in f.what
               for f in found)


def test_check_reports_a_clone_covering_one_chain() -> None:
    rec = derive(derive(records(), PMHC), EPITOPE)
    ch = (derive(chains(), CLONOTYPE)
          .with_columns(pl.lit("CXorphan").alias(CLONE.column)).drop("species"))
    found = checks.check({"records": rec, "chains": ch})
    assert any("two chains of different gene" in f.what for f in found)


# --------------------------------------------------------------------------------------------
# Lifecycle
# --------------------------------------------------------------------------------------------

def _present(*ids: tuple[str, str]) -> pl.DataFrame:
    return pl.DataFrame({"id": [i for i, _ in ids], "level": [lv for _, lv in ids]},
                        schema={"id": pl.Utf8, "level": pl.Utf8})


def test_a_first_release_names_itself_as_both_endpoints() -> None:
    out = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    row = out.row(0, named=True)
    assert row["state"] == lifecycle.ACTIVE
    assert row["first_release"] == row["last_release"] == "v1"


def test_an_id_that_survives_keeps_its_first_release_and_moves_its_last() -> None:
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    v2 = lifecycle.advance(v1, _present(("CTaaa", "clonotype")), release="v2")
    row = v2.row(0, named=True)
    assert row["first_release"] == "v1"
    assert row["last_release"] == "v2"


def test_a_retirement_freezes_the_release_that_last_carried_the_id() -> None:
    """The question a consumer asks is which release still had it, so `last_release` must not move."""
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    v2 = lifecycle.advance(v1, _present(("CTbbb", "clonotype")), release="v2")
    gone = v2.filter(pl.col("id") == "CTaaa").row(0, named=True)
    assert gone["state"] == lifecycle.RETIRED
    assert gone["last_release"] == "v1"
    v3 = lifecycle.advance(v2, _present(("CTbbb", "clonotype")), release="v3")
    assert v3.filter(pl.col("id") == "CTaaa").row(0, named=True)["last_release"] == "v1"


def test_a_returning_id_is_active_again_with_its_original_first_release() -> None:
    """A derived id *is* the hash of its key, so the same key returning has to give the same id."""
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    v2 = lifecycle.advance(v1, _present(("CTbbb", "clonotype")), release="v2")
    v3 = lifecycle.advance(v2, _present(("CTaaa", "clonotype")), release="v3")
    back = v3.filter(pl.col("id") == "CTaaa").row(0, named=True)
    assert back["state"] == lifecycle.ACTIVE
    assert back["first_release"] == "v1"
    assert back["last_release"] == "v3"
    report = lifecycle.compare(v2, _present(("CTaaa", "clonotype")), release="v3")
    assert report.returned.height == 1
    assert report.added.is_empty()


def test_a_replacement_is_recorded_only_where_the_caller_names_one() -> None:
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    v2 = lifecycle.advance(v1, _present(("CTbbb", "clonotype")), release="v2",
                           replaced_by={"CTaaa": "CTbbb"})
    assert v2.filter(pl.col("id") == "CTaaa").row(0, named=True)["replaced_by"] == "CTbbb"
    assert v2.filter(pl.col("id") == "CTbbb").row(0, named=True)["replaced_by"] == ""


def test_the_lifecycle_survives_a_round_trip_through_a_file(tmp_path) -> None:
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype"),
                                                       ("PMbbb", "pmhc")), release="v1")
    path = tmp_path / "lifecycle.tsv"
    lifecycle.write(v1, path)
    assert lifecycle.read(path).equals(v1)
    assert lifecycle.read(tmp_path / "absent.tsv").is_empty()


def test_resolve_finds_an_id_and_says_nothing_about_one_it_never_saw() -> None:
    v1 = lifecycle.advance(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    assert lifecycle.resolve(v1, "CTaaa")["level"] == "clonotype"
    assert lifecycle.resolve(v1, "CTzzz") is None


def test_an_unknown_prefix_is_reported() -> None:
    assert lifecycle.unknown_prefixes(_present(("ZZaaa", "clonotype"))) == ["ZZaaa"]
    assert lifecycle.unknown_prefixes(_present(("CTaaa", "clonotype"))) == []


def test_the_counts_table_names_every_level_even_at_zero() -> None:
    report = lifecycle.compare(lifecycle.empty(), _present(("CTaaa", "clonotype")), release="v1")
    counts = report.counts()
    assert counts.height == len(LEVELS)
    assert counts.filter(pl.col("level") == "clonotype")["added"][0] == 1
    assert counts.filter(pl.col("level") == "pmhc")["added"][0] == 0
