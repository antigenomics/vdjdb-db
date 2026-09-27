"""The legacy k-mer V/J guesser, which is still in the build path.

`arda` repairs a CDR3 against a named segment; it does not propose one. Chunks may leave `v.alpha`
or `j.beta` blank, so `_legacy_guess.guess_segments` fills them from the germline parts in `res/`.
Dropping this guess cost 13,844 of 284,546 rows of `vdjdb.txt` when it was tried, which is why the
scan survives until the Pgen guesser of #462.

It had no test, while three modules call it: `_legacy_guess`, `assemble/master` and
`curate/nomenclature`.

Note on measuring it: `_legacy_fixer/__init__.py` puts its own directory on `sys.path` and imports
`Cdr3Fixer` as a top-level module, the way `py_src/` did. `--cov=vdjdb` filters by import name, so
it reports these files at 0 % however much they run. Measure by path instead:

    uv run pytest tests/unit/test_legacy_guess.py --cov=src/vdjdb/annotate/_legacy_fixer

which reads 69 % with these tests (Cdr3Fixer.py 74 %).
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb._legacy_guess import guess_segments
from vdjdb.config import Paths

pytestmark = pytest.mark.skipif(
    not (Paths.discover().res / "segments.txt").exists(),
    reason="res/segments.txt is required by the legacy guesser")


def _keys(rows):
    return pl.DataFrame(rows, schema={"species": pl.Utf8, "cdr3": pl.Utf8,
                                      "v": pl.Utf8, "j": pl.Utf8}, orient="row")


@pytest.fixture(scope="module")
def fixer():
    from vdjdb.annotate._legacy_fixer import Cdr3Fixer
    res = Paths.discover().res
    return Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))


def test_a_given_segment_is_never_replaced():
    """`v` and `j` are the join key back to the table; overwriting them fans the output out."""
    out = guess_segments(_keys([("HomoSapiens", "CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-7*01")]))
    assert out["__gv"][0] == "TRBV19*01"
    assert out["__gj"][0] == "TRBJ2-7*01"
    assert out["v"][0] == "TRBV19*01" and out["j"][0] == "TRBJ2-7*01"


@pytest.mark.xfail(strict=True, reason=(
    "Cdr3Fixer.guess_id puts `return \"\"` inside the five-prime loop, so it tries one prefix and "
    "gives up; the Groovy it was ported from (attic/Cdr3Fixer.groovy:165-171) returns after the "
    "loop. Measured 2026-09-27: 711 chunk rows carry a CDR3 with no V call, over 644 distinct "
    "(species, cdr3) keys, and the corrected loop names a V for 640 of them. Fixing it changes "
    "build output, so it needs an entry in rules/expected_diffs.toml first."))
def test_a_blank_v_is_guessed_from_the_cdr3():
    out = guess_segments(_keys([("HomoSapiens", "CASSIRSSYEQYF", "", "TRBJ2-7*01")]))
    assert out["__gv"][0].startswith("TRBV"), out["__gv"][0]
    assert out["v"][0] == "", "the join key must stay blank"


def test_a_blank_v_currently_stays_blank_without_raising():
    """What the five-prime bug above means in practice: no call, and no exception."""
    out = guess_segments(_keys([("HomoSapiens", "CASSIRSSYEQYF", "", "TRBJ2-7*01")]))
    assert out["__gv"][0] == ""
    assert out["v"][0] == ""


def test_a_blank_j_is_guessed_from_the_cdr3():
    out = guess_segments(_keys([("HomoSapiens", "CASSIRSSYEQYF", "TRBV19*01", "")]))
    assert out["__gj"][0].startswith("TRBJ"), out["__gj"][0]


def test_the_chain_is_read_off_whichever_segment_is_present():
    """With no `gene` given, the locus comes from the segment that is there."""
    alpha = guess_segments(_keys([("HomoSapiens", "CAVRDSNYQLIW", "TRAV3*01", "")]))
    assert alpha["__gj"][0].startswith("TRAJ"), alpha["__gj"][0]


def test_an_explicit_gene_overrides_the_sniffed_one():
    """With both segments blank there is nothing to sniff, so `gene` is the only signal."""
    out = guess_segments(_keys([("HomoSapiens", "CAVRDSNYQLIW", "", "")]), gene="alpha")
    assert out["__gj"][0].startswith("TRAJ"), out["__gj"][0]


def test_nothing_to_guess_short_circuits_without_loading_the_reference():
    """No blank means no `Cdr3Fixer`, which is what keeps a full build from paying for it."""
    from vdjdb.annotate.cdr3fix import guess_missing_segments
    keys = _keys([("HomoSapiens", "CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-7*01")])
    out = guess_missing_segments(keys)
    assert out["__gv"].to_list() == ["TRBV19*01"]
    assert out["__gj"].to_list() == ["TRBJ2-7*01"]


def test_an_unguessable_cdr3_yields_a_blank_not_an_exception():
    out = guess_segments(_keys([("HomoSapiens", "XXXXXXXXXX", "", "")]), gene="beta")
    assert out["__gv"][0] == "" and out["__gj"][0] == ""


def test_the_species_key_is_case_sensitive():
    """`guess_id` looks up `<species>.<gene>`, and the parts table spells it `HomoSapiens`.

    `_load_segments_sequence_data` does not lowercase where `_load_segments_data` does, so a caller
    passing `homosapiens` gets no guess at all. VDJdb's own spelling is the CamelCase one, which is
    why the build path works; anything reaching for this directly must match it.
    """
    from vdjdb.annotate._legacy_fixer import Cdr3Fixer
    res = Paths.discover().res
    fx = Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))
    assert fx.guess_id("CASSIRSSYEQYF", "HomoSapiens", "beta", False) == "TRBJ2-7"
    assert fx.guess_id("CASSIRSSYEQYF", "homosapiens", "beta", False) == ""


def test_get_closest_id_resolves_an_allele_free_call(fixer):
    """Chunks routinely write `TRBV19` with no allele; the fixer picks a concrete one."""
    got = fixer.get_closest_id("homosapiens", "TRBV19")
    assert got.startswith("TRBV19"), got


def test_get_segment_seq_returns_none_for_an_unknown_segment(fixer):
    assert fixer.get_segment_seq("homosapiens", "TRBV999*01") is None


def test_fix_reports_a_bad_segment_rather_than_raising(fixer):
    """An unknown segment is a reported fix type, not an exception."""
    result = fixer.fix("CASSIRSSYEQYF", "TRBV999*01", "homosapiens", True)
    assert result.FixType.name == "FailedBadSegment"
    assert result.cdr3 == "CASSIRSSYEQYF", "the sequence is returned unchanged"
