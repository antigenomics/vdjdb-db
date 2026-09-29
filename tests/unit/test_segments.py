"""The V/J proposal (#462, #658): filled only where the record named nothing, and never overwriting it.

`propose` answers one question - what segment does this junction name, given that its record does not -
and both readers of the answer go through it: `annotate.cdr3fix.markup` hands the J to the repair, and
`chains` reports both as `v.inferred` / `j.inferred`.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.annotate import segments
from vdjdb.annotate.cdr3fix import ensure_reference


def _keys(rows: list[tuple[str, str, str, str]]) -> pl.DataFrame:
    return pl.DataFrame(rows, schema=["species", "cdr3", "v", "j"], orient="row")


@pytest.fixture(scope="module")
def proposed() -> pl.DataFrame:
    ensure_reference()
    return segments.propose(_keys([
        ("HomoSapiens", "CASSIRSSYEQYF", "TRBV10-3*01", "TRBJ2-7*01"),   # nothing to do
        ("HomoSapiens", "CASSIRSSYEQYF", "", "TRBJ2-7*01"),              # V missing
        ("HomoSapiens", "CASSIRSSYEQYF", "TRBV10-3*01", ""),             # J missing
        ("HomoSapiens", "CASSIRSSYEQYF", "", ""),                        # both missing
        ("RattusNorvegicus", "CASSIRSSYEQYF", "", ""),                   # neither source
    ]), "beta").sort("species", "v", "j")


def _row(frame: pl.DataFrame, v: str, j: str, species: str = "HomoSapiens") -> dict:
    return frame.filter((pl.col("species") == species) & (pl.col("v") == v)
                        & (pl.col("j") == j)).row(0, named=True)


def test_a_call_the_record_carries_is_never_replaced(proposed):
    """`v` and `j` are the join key back to the table, so the proposal goes in its own columns."""
    r = _row(proposed, "TRBV10-3*01", "TRBJ2-7*01")
    assert r["__gv"] == "TRBV10-3*01" and r["__gj"] == "TRBJ2-7*01"


def test_the_missing_side_is_the_only_one_proposed(proposed):
    v_missing = _row(proposed, "", "TRBJ2-7*01")
    assert v_missing["__gv"].startswith("TRBV") and v_missing["__gj"] == "TRBJ2-7*01"
    j_missing = _row(proposed, "TRBV10-3*01", "")
    assert j_missing["__gj"].startswith("TRBJ") and j_missing["__gv"] == "TRBV10-3*01"


def test_both_sides_are_proposed_when_neither_is_named(proposed):
    r = _row(proposed, "", "")
    assert r["__gv"].startswith("TRBV") and r["__gj"].startswith("TRBJ")


def test_a_species_with_no_model_is_answered_from_ardas_germline_table():
    """The second source, and the only one macaque has: there is no macaque recombination model, and
    the retired k-mer scanner named nothing for its 74 blank J calls either. arda carries a rhesus
    anchor table - 73 TRAV, 69 TRAJ, 117 TRBV, 17 TRBJ alleles - so the germline can answer.

    Alpha, because all 74 of those blanks are `j.alpha`.
    """
    ensure_reference()
    r = segments.propose(_keys([("MacacaMulatta", "CAVSGGYQKVTF", "", "")]), "alpha").row(0, named=True)
    assert r["__gj"].startswith("TRAJ"), r
    assert r["__gv"].startswith("TRAV"), r


def test_a_species_with_neither_source_gets_nothing_rather_than_another_organisms_answer(proposed):
    """arda ships no rat TR anchors and vdjtools no rat model, so a rat chain is left blank. Naming a
    human segment on it would be worse than naming none.
    """
    r = _row(proposed, "", "", species="RattusNorvegicus")
    assert r["__gv"] == "" and r["__gj"] == ""


def test_no_key_is_lost_or_duplicated(proposed):
    assert proposed.height == 5
    assert proposed.select("species", "cdr3", "v", "j").is_duplicated().sum() == 0


def test_a_frame_with_nothing_to_propose_is_returned_unchanged():
    keys = _keys([("HomoSapiens", "CASSIRSSYEQYF", "TRBV10-3*01", "TRBJ2-7*01")])
    out = segments.propose(keys, "beta")
    assert out["__gv"].to_list() == ["TRBV10-3*01"] and out["__gj"].to_list() == ["TRBJ2-7*01"]


def test_the_locus_is_read_off_whichever_call_is_present_when_the_caller_does_not_say():
    """A frame mixing both loci is what `chains` hands over; `markup` always knows its chain word."""
    ensure_reference()
    mixed = _keys([("HomoSapiens", "CAVRDSNYQLIW", "TRAV3*01", ""),
                   ("HomoSapiens", "CASSIRSSYEQYF", "TRBV10-3*01", "")])
    out = segments.propose(mixed).sort("v")
    assert out.filter(pl.col("v") == "TRAV3*01")["__gj"].item().startswith("TRAJ")
    assert out.filter(pl.col("v") == "TRBV10-3*01")["__gj"].item().startswith("TRBJ")


def test_the_germline_call_needs_three_residues_of_agreement():
    """Two residues of a junction's end match most J alleles of a locus, so a shorter run would name
    whichever allele sorted first rather than whichever fits.
    """
    anchors = (("TRBJ2-7*01", "SYEQYF"), ("TRBJ1-1*01", "NTEAFF"))
    assert segments._germline_call("CASSIRSSYEQYF", anchors, five_prime=False) == "TRBJ2-7*01"
    assert segments._germline_call("CASSIRSSYQQQQ", anchors, five_prime=False) == ""


def test_a_functional_allele_wins_a_tie_with_a_pseudogene():
    """A pseudogene cannot be the segment of an expressed receptor, so where two alleles agree with
    the junction equally far the functional one is the answer. `_anchors` encodes this in its sort,
    which is also what makes the tiebreak the same on every host (rule 7).
    """
    ensure_reference()
    from arda.cdr3fix import load_anchors

    table = load_anchors("human")
    for locus, segment in (("TRB", "V"), ("TRB", "J"), ("TRA", "V"), ("TRA", "J")):
        rows = segments._anchors("human", locus, segment)
        assert rows, f"no {locus}{segment} anchors"
        functional = [table[(segment, allele)].functionality == "F" for allele, _t in rows]
        assert functional == sorted(functional, reverse=True), f"{locus}{segment} is not ordered"
