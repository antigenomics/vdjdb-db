"""Coordinate conversions: every pair round-trips, and junction is never confused with CDR3."""
from __future__ import annotations

import polars as pl
import pytest
from hypothesis import given
from hypothesis import strategies as st

from vdjdb.convert import coords

AA = st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=0, max_size=30)


@given(AA)
def test_junction_and_cdr3_differ_by_exactly_the_two_anchors(cdr3):
    """The whole reason this module exists: `junction_aa` is two residues longer than `cdr3_aa`."""
    junction = coords.junction_from_cdr3(cdr3)
    if not cdr3:
        assert junction == ""
        return
    assert len(junction) == len(cdr3) + 2
    assert coords.cdr3_from_junction(junction) == cdr3


@pytest.mark.parametrize("junction,expected", [
    ("CASSIRSSYEQYF", "ASSIRSSYEQY"),
    ("CASSIRSSYEQYW", "ASSIRSSYEQY"),   # Trp118, not only Phe
    ("CAF", "A"),                       # the shortest junction that has a CDR3 at all
    ("CF", ""), ("C", ""), ("", ""),    # anchors only, or less: no CDR3, not a negative slice
])
def test_stripping_the_anchors(junction, expected):
    assert coords.cdr3_from_junction(junction) == expected


@given(AA)
def test_the_polars_expression_agrees_with_the_scalar_function(junction):
    got = pl.DataFrame({"cdr3": [junction]}).select(coords.cdr3_aa()).item()
    assert got == coords.cdr3_from_junction(junction)


@given(st.integers(min_value=0, max_value=10_000))
def test_origin_round_trips(i):
    assert coords.to_zero_based(coords.to_one_based(i)) == i


@given(st.integers(min_value=1, max_value=10_000))
def test_closedness_round_trips(end):
    assert coords.close_end(coords.open_end(end)) == end
    assert coords.open_end(coords.close_end(end)) == end


@given(st.integers(min_value=0, max_value=3_000))
def test_amino_acid_to_nucleotide_round_trips_zero_based(i):
    assert coords.nt_to_aa(coords.aa_to_nt(i)) == i


@given(st.integers(min_value=1, max_value=3_000))
def test_amino_acid_to_nucleotide_round_trips_one_based(i):
    assert coords.nt_to_aa(coords.aa_to_nt(i, one_based=True), one_based=True) == i


def test_the_two_origins_do_not_agree_by_accident():
    """If these ever coincided, a missing `one_based=` would be undetectable."""
    assert coords.aa_to_nt(1) == 3 and coords.aa_to_nt(1, one_based=True) == 1


@given(st.integers(min_value=0, max_value=500), st.integers(min_value=0, max_value=5_000))
def test_junction_space_to_sequence_space_round_trips(pos, start):
    assert coords.from_sequence(coords.to_sequence(pos, start), start) == pos


# -- the fitted V/J boundary conversion ---------------------------------------------------------

@pytest.mark.parametrize("nt", range(0, 40))
def test_the_boundary_conversion_has_one_definition(nt: int) -> None:
    """The scalar and the column form are the same arithmetic, asserted rather than assumed.

    Two spellings of one coordinate conversion is what the module comment warns about, and this is
    the only thing that stops them drifting - the scalar is what a reader reasons with and the
    expression is what the build runs on 187,055 boundaries.
    """
    scalar = coords.nt_to_aa_boundary(nt)
    column = (pl.DataFrame({"nt": [nt]})
              .select(coords.nt_to_aa_boundary_expr(pl.col("nt")))
              .item())
    assert scalar == column


def test_the_boundary_conversion_ceilings_where_nt_to_aa_floors() -> None:
    """Fitted, not reasoned: each is the count of residues the segment touches, so a partly-covered
    codon counts. Under floor the exact-match rates against the external truth set in
    `tests/release/test_cdr3fix_accuracy.py` drop by tens of percent.
    """
    assert [coords.nt_to_aa_boundary(i) for i in range(7)] == [0, 1, 1, 1, 2, 2, 2]
    assert [coords.nt_to_aa(i) for i in range(7)] == [0, 0, 0, 1, 1, 1, 2]


@pytest.mark.parametrize("bad", [None, -1, -7])
def test_an_unmapped_boundary_stays_unmapped(bad) -> None:
    """A null or a negative is "the engine declined", and it must not become residue 0."""
    got = (pl.DataFrame({"nt": [bad]}, schema={"nt": pl.Int64})
           .select(coords.nt_to_aa_boundary_expr(pl.col("nt"))).item())
    assert got == -1
