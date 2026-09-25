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


# -- allele disambiguation from the CDR3 (#327) --------------------------------------------------

def _traj24(calls, cdr3s, species=None):
    n = len(calls)
    return pl.DataFrame({
        "species": species or ["HomoSapiens"] * n,
        "j.alpha": calls, "cdr3.alpha": cdr3s,
        "v.alpha": [""] * n, "v.beta": [""] * n, "d.beta": [""] * n, "j.beta": [""] * n,
    })


def test_the_cdr3_decides_the_allele_whatever_the_submitter_wrote():
    """73 of the 111 records explicitly called *01 carry the *02 signature; none carry *01's."""
    df = _traj24(["TRAJ24", "TRAJ24*01", "TRAJ24*02"], ["CAWGKLQF"] * 3)
    out, report = N.disambiguate_alleles(df)
    assert out["j.alpha"].to_list() == ["TRAJ24*02"] * 3
    assert report["rows"].sum() == 2          # the one already correct is not a change


def test_the_other_allele_signature_is_honoured_symmetrically():
    """`WGKFEF` appears zero times in the corpus today; the rule must still be the right one."""
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24*02"], ["CAWGKFEF"]))
    assert out["j.alpha"].to_list() == ["TRAJ24*01"]


def test_a_cdr3_with_no_signature_is_left_alone():
    """364 records have a CDR3 trimmed short of the anchor. No evidence, no correction."""
    out, report = N.disambiguate_alleles(_traj24(["TRAJ24", "TRAJ24*01"], ["CAVSDLE", "CAVSDLE"]))
    assert out["j.alpha"].to_list() == ["TRAJ24", "TRAJ24*01"]
    assert report.is_empty()


def test_a_cdr3_carrying_both_signatures_contradicts_itself_and_is_refused():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24"], ["CAWGKLQFWGKFEF"]))
    assert out["j.alpha"].to_list() == ["TRAJ24"]


def test_another_species_is_not_touched_by_a_human_rule():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24"], ["CAWGKLQF"], ["MusMusculus"]))
    assert out["j.alpha"].to_list() == ["TRAJ24"]


def test_a_different_gene_is_not_touched():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ42*01"], ["CAWGKLQF"]))
    assert out["j.alpha"].to_list() == ["TRAJ42*01"]


def test_the_allele_rename_carries_its_evidence_into_the_ledger():
    """Not injective on value alone: the fixer resolves a bare `TRAJ24` to `*01`, so the reference
    ships the same value for the records the CDR3 corrects and the ones it does not."""
    report = pl.DataFrame({"issue": ["#327"] * 2, "column": ["j.alpha"] * 2,
                           "species": ["HomoSapiens"] * 2, "from": ["TRAJ24", "TRAJ24*01"],
                           "to": ["TRAJ24*02"] * 2, "signature": ["WGKLQF"] * 2,
                           "rows": pl.Series([974, 73], dtype=pl.UInt32)})
    block = N.render_allele_renames(report, lambda sp, c: "TRAJ24*01" if c == "TRAJ24" else c)
    assert block.count("[[rename]]") == 1, "both rows resolve to one reference value"
    assert 'from = "TRAJ24*01"' in block and 'to = "TRAJ24*02"' in block
    assert 'when_contains = "WGKLQF"' in block
    assert "records = 1047" in block
