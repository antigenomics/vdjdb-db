"""The V/J guesser (#462): filled only where curation left a gap, and never overwriting it."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.annotate import segments

RECORDS = pl.DataFrame({"record_id": ["r1", "r2", "r3", "r4"],
                        "species": ["HomoSapiens"] * 3 + ["MacacaMulatta"]})
CHAINS = pl.DataFrame({
    "record_id": ["r1", "r2", "r3", "r4"],
    "gene": ["TRB"] * 4,
    "cdr3": ["CASSIRSSYEQYF"] * 4,
    "v.segm": ["TRBV10-3*01", "", "TRBV10-3*01", ""],
    "j.segm": ["TRBJ2-7*01", "TRBJ2-7*01", "", ""],
})


@pytest.fixture(scope="module")
def guessed() -> pl.DataFrame:
    return segments.add_inferred_segments(CHAINS, RECORDS)


def test_a_curated_call_is_never_overwritten(guessed):
    """The curated columns are the publication's. This stage only fills a gap beside them."""
    assert guessed["v.segm"].to_list() == CHAINS["v.segm"].to_list()
    assert guessed["j.segm"].to_list() == CHAINS["j.segm"].to_list()
    r1 = guessed.filter(pl.col("record_id") == "r1").row(0, named=True)
    assert r1["v.inferred"] == "" and r1["j.inferred"] == ""


def test_the_missing_side_is_the_only_one_filled(guessed):
    r2 = guessed.filter(pl.col("record_id") == "r2").row(0, named=True)
    assert r2["v.inferred"].startswith("TRBV") and r2["j.inferred"] == ""
    r3 = guessed.filter(pl.col("record_id") == "r3").row(0, named=True)
    assert r3["j.inferred"].startswith("TRBJ") and r3["v.inferred"] == ""


def test_a_species_with_no_model_gets_nothing_rather_than_a_wrong_organisms_answer(guessed):
    r4 = guessed.filter(pl.col("record_id") == "r4").row(0, named=True)
    assert r4["v.inferred"] == "" and r4["j.inferred"] == ""


def test_no_row_is_dropped(guessed):
    assert guessed.height == CHAINS.height
    assert "species" not in guessed.columns


def test_the_legacy_kmer_v_guesser_is_the_thing_this_replaces():
    """`Cdr3Fixer.guess_id` puts `return ""` inside the five-prime loop, so it tries one prefix
    length and gives up: 3 non-empty V guesses in 4,000 sequences against 3,797 for J.

    Asserted rather than described, because it is the whole justification for #462 and because a
    future reader will otherwise assume the k-mer scan merely performs worse.
    """
    from vdjdb.annotate._legacy_fixer import Cdr3Fixer
    from vdjdb.config import Paths

    res = Paths.discover().res
    fx = Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))
    seqs = ["CASSIRSSYEQYF", "CASSLAPGATNEKLFF", "CASSPGQGAYEQYF", "CASSQDRGNTGELFF",
            "CASSEGWHSYEQYF", "CASSLGQAYEQYF"]
    v = [fx.guess_id(s, "HomoSapiens", "beta", True) for s in seqs]
    j = [fx.guess_id(s, "HomoSapiens", "beta", False) for s in seqs]
    assert not any(v), "the five-prime branch returns early; if this passes guesses, it was fixed"
    assert all(j), "the three-prime branch works, which is what makes the V result a bug"
