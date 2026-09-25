"""IMGT segment nomenclature (#389): mechanical respellings only, never a guess."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate import nomenclature as N


@pytest.mark.parametrize("call,expected", [
    # the /DV genes IMGT names once for both loci -- 1,377 chains write the short form
    ("TRAV14", "TRAV14/DV4"),
    ("TRAV29", "TRAV29/DV5"),
    ("TRAV29DV5", "TRAV29/DV5"),        # the slash dropped entirely
    ("TRBJ1.2", "TRBJ1-2"),             # a dot where a dash belongs
    ("TRBJ 2-7", "TRBJ2-7"),            # a space
    ("TRBD1-1*01", "TRBD1*01"),         # a D gene with a suffix IMGT does not use
    ("TCRBD2*02", "TRBD2*02"),          # the pre-IMGT TCR prefix
])
def test_mechanical_respellings(call, expected):
    assert N.normalise_call(call, "HomoSapiens") == expected


@pytest.mark.parametrize("call", [
    "TRBV10-3*01",      # already IMGT
    "TRBV8",            # several candidates: TRBV8-1, TRBV8-2 -- a curation question, not a fix
    "TRAJ16.5",         # no IMGT counterpart at all
    "",
])
def test_a_call_that_is_correct_or_undecidable_is_left_alone(call):
    assert N.normalise_call(call, "HomoSapiens") is None


def test_a_species_with_no_authority_is_never_touched():
    """An unchecked rewrite is worse than an odd spelling."""
    assert N.normalise_call("TRAV14", "GallusGallus") is None


def test_a_multi_call_is_normalised_member_by_member_and_sorted():
    """`TRBD2,TRBD1` and `TRBD1,TRBD2` are two spellings of one fact -- 89 chains' worth."""
    assert N.normalise_call("TRBD2,TRBD1", "HomoSapiens") == "TRBD1,TRBD2"
    assert N.normalise_call("TRBD2*01 or TCRBD2*02", "HomoSapiens") == "TRBD2*01,TRBD2*02"


def test_a_multi_call_with_one_unresolvable_member_is_refused_whole():
    """Half a correction is not a correction."""
    assert N.normalise_call("TRBD1,TRBV8", "HomoSapiens") is None


def test_the_frame_is_rewritten_and_the_report_accounts_for_every_change():
    df = pl.DataFrame({
        "species": ["HomoSapiens", "HomoSapiens", "MusMusculus"],
        "v.alpha": ["TRAV14", "TRAV14", "TRAV21-DV12"],
        "j.alpha": ["TRAJ3*01", "TRAJ3*01", "TRAJ3*01"],
        "v.beta": ["", "", ""], "d.beta": ["", "", ""], "j.beta": ["", "", ""],
    })
    out, report = N.harmonise_segments(df)
    assert out["v.alpha"].to_list() == ["TRAV14/DV4", "TRAV14/DV4", "TRAV21/DV12"]
    assert report["rows"].sum() == 3
    assert set(report["from"]) == {"TRAV14", "TRAV21-DV12"}
    assert report.filter(pl.col("from") == "TRAV14")["rows"][0] == 2


def test_nothing_to_do_yields_an_empty_report_not_an_error():
    df = pl.DataFrame({"species": ["HomoSapiens"], "v.alpha": ["TRAV12-2*01"],
                       "j.alpha": [""], "v.beta": [""], "d.beta": [""], "j.beta": [""]})
    out, report = N.harmonise_segments(df)
    assert report.is_empty() and out["v.alpha"][0] == "TRAV12-2*01"


# -- the generated ledger block ----------------------------------------------------------------

def test_a_rename_is_declared_only_when_the_fixer_leaves_the_old_spelling_alone():
    """The injectivity condition, and the reason the ledger can trust the block.

    `get_closest_id` simplifies `TRAV6-7-DV9` to `TRAV6` and then tries `TRAV6-1*01`, `TRAV6-2*01`,
    ... taking the first hit -- so the reference ships `TRAV6-1*01` and is indistinguishable from
    records that genuinely are TRAV6-1. Declaring that rename rewrites both.
    """
    report = pl.DataFrame({"column": ["v.alpha", "v.alpha"], "species": ["MusMusculus"] * 2,
                           "from": ["TRAV6-7-DV9", "TRAV14"], "to": ["TRAV6-7/DV9", "TRAV14/DV4"],
                           "rows": pl.Series([15, 2], dtype=pl.UInt32)})
    resolve = lambda sp, call: {                     # noqa: E731
        "TRAV6-7-DV9": "TRAV6-1*01", "TRAV6-7/DV9": "TRAV6-7/DV9*01",
        "TRAV14": "TRAV14", "TRAV14/DV4": "TRAV14/DV4*01"}.get(call, call)
    block = N.render_renames(report, resolve)
    assert "TRAV6-7-DV9" not in block and "TRAV6-1" not in block
    assert 'from = "TRAV14"' in block and 'to = "TRAV14/DV4*01"' in block


def test_a_multi_call_is_never_declared_as_a_rename():
    """`fix_both` splits it and keeps the best member, so what ships is a selection, not a rename."""
    report = pl.DataFrame({"column": ["d.beta"], "species": ["HomoSapiens"],
                           "from": ["TRBD2,TRBD1"], "to": ["TRBD1,TRBD2"],
                           "rows": pl.Series([43], dtype=pl.UInt32)})
    assert "[[rename]]" not in N.render_renames(report, lambda sp, c: c)


def test_the_generated_block_is_replaced_in_place_not_appended(tmp_path):
    path = tmp_path / "rules.toml"
    path.write_text('[[rule]]\nid = "keep-me"\n')
    report = pl.DataFrame({"column": ["v.alpha"], "species": ["HomoSapiens"],
                           "from": ["TRAV14"], "to": ["TRAV14/DV4"],
                           "rows": pl.Series([2], dtype=pl.UInt32)})
    N.write_renames(report, path, lambda sp, c: c)
    N.write_renames(report, path, lambda sp, c: c)
    text = path.read_text()
    assert text.count("[[rename]]") == 1 and "keep-me" in text
