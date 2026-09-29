"""The junction-anchor check: what it flags, what it proposes, and what it refuses to call a defect."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate import anchors


def _master(rows: list[dict]) -> pl.DataFrame:
    """A master-shaped frame with only the columns the check reads."""
    base = {"record_id": "", "species": "HomoSapiens", "chunk.file": "PMID_1.txt", "chunk.row": 1,
            "cdr3.alpha": "", "v.alpha": "", "j.alpha": "", "__cdr3old.alpha": None,
            "__vcanon.alpha": True, "__jcanon.alpha": True,
            "cdr3.beta": "", "v.beta": "", "j.beta": "", "__cdr3old.beta": None,
            "__vcanon.beta": True, "__jcanon.beta": True}
    return pl.DataFrame([{**base, **r} for r in rows])


def test_the_germline_anchor_is_read_per_segment_and_not_assumed_to_be_f_or_w():
    """`TRAJ35*01` is `IGFGNVLHC` and mouse `TRAJ47*01` is `HYANKMIC`.

    A hardcoded "ends with F or W" test calls a junction on either of those broken. Measured, that
    mistake flags 125 correct rows.
    """
    assert anchors.templated("HomoSapiens", "J", "TRAJ35*01") == "IGFGNVLHC"
    assert anchors.templated("MusMusculus", "J", "TRAJ47*01") == "HYANKMIC"
    assert (anchors.templated("HomoSapiens", "V", "TRBV6-1*01") or "").startswith("C")
    # An allele-less call falls back to *01 rather than giving up: a record naming no allele is not
    # claiming one.
    assert anchors.templated("HomoSapiens", "J", "TRBJ1-4") == "TNEKLFF"


def test_a_species_or_segment_the_reference_does_not_have_is_unchecked_not_broken():
    assert anchors.templated("Nonexistent", "V", "TRBV6-1*01") is None
    assert anchors.templated("HomoSapiens", "V", "TRBV999-9*01") is None
    assert anchors.templated("HomoSapiens", "V", "") is None


def test_a_canonical_junction_is_not_flagged():
    defect, repair = anchors.classify("CASSNEKLFF", "HomoSapiens", "TRBV6-1*01", "TRBJ1-4*01")
    assert (defect, repair) == ("ok", None)


def test_a_missing_j_anchor_is_named_and_the_germline_residue_is_proposed():
    """The reported case: `TRBJ1-4` is `TNEKLFF`, so a junction ending `NEKLF` is one Phe short."""
    defect, repair = anchors.classify("CASSNEKLF", "HomoSapiens", "TRBV6-1*01", "TRBJ1-4*01")
    assert defect == "J absent anchor"
    assert repair == "CASSNEKLFF"


def test_both_anchors_missing_are_repaired_together():
    """`ASSNEKLF` is short a Cys in front and a Phe behind, and the two repairs must compose."""
    defect, repair = anchors.classify("ASSNEKLF", "HomoSapiens", "TRBV6-1*01", "TRBJ1-4*01")
    assert defect == "V absent anchor, J absent anchor"
    assert repair == "CASSNEKLFF"


def test_a_mis_read_cys104_is_substituted_rather_than_prepended():
    """No TCR folds without Cys104, so a first residue that is not Cys where the body aligns is a
    read error, not a variant - and prepending would leave the wrong residue in place."""
    for wrong in ("GASSNEKLFF", "WASSNEKLFF", "FASSNEKLFF"):
        defect, repair = anchors.classify(wrong, "HomoSapiens", "TRBV6-1*01", "TRBJ1-4*01")
        assert defect == "V corrupt anchor", wrong
        assert repair == "CASSNEKLFF", wrong


def test_a_short_germline_contribution_does_not_stop_the_trim_being_recognised():
    """A V contributes three to five junction residues before the N region takes over - `TRAV27*01`
    is `CAG` - so an under-trim cannot be recognised by deep germline agreement. `YLCSSQEGGYGYTFGSG`
    on `TRBV29-1*01` (`CSVE`) aligns 0 residues as given and 2 from its Cys, and 2 > 0 is the whole
    signal."""
    defect, repair = anchors.classify("YLCSSQEGGYGYTFGSG", "HomoSapiens",
                                      "TRBV29-1*01", "TRBJ1-2*01")
    assert defect == "V under-trimmed, J under-trimmed"
    assert repair == "CSSQEGGYGYTF"


def test_framework_carried_past_an_anchor_is_trimmed_at_both_ends():
    """`YFC...` in front and `...FGXG` behind: the junction plus the V and J framework around it.

    `arda.cdr3fix` trims one end and leaves the other, which is why these reach the report at all.
    """
    defect, repair = anchors.classify("YFCASSNEKLFFGSG", "HomoSapiens", "TRBV6-1*01", "TRBJ1-4*01")
    assert defect == "V under-trimmed, J under-trimmed"
    assert repair == "CASSNEKLFF"


def test_a_junction_ending_in_the_unusual_germline_residue_is_correct():
    """`TRAJ35*01` encodes Cys at 118, so `...GFGNVLHC` is the canonical junction, not a defect."""
    defect, repair = anchors.classify("CAASLGFGNVLHC", "HomoSapiens", "TRAV13-1*01", "TRAJ35*01")
    assert defect == "ok"
    assert repair is None


def test_a_universal_anchor_against_a_disagreeing_table_blames_the_table():
    """Mouse `TRAJ47*01` is `HYANKMIC` and every record reads `DYANKMIF`; `YANKMI` is identical, so
    the frame is right and both terminal residues are not. 95 rows, and the evidence points at the
    reference this repository does not own - so it is reported apart and never repaired."""
    defect, repair = anchors.classify("CPDYANKMIF", "MusMusculus", "TRAV8D-1*01", "TRAJ47*01")
    assert defect == "J anchor table suspect"
    assert repair is None


def test_a_defect_the_germline_does_not_explain_gets_no_repair():
    """Saying "broken, and here is a guess" would be worse than saying "broken"."""
    defect, repair = anchors.classify("CASSPLPGT", "HomoSapiens", "TRBV11-2*01", "TRBJ2-1*01")
    assert defect == "J unexplained"
    assert repair is None


def test_the_repair_is_proposed_against_the_submitted_sequence_not_the_shipped_one():
    """An edit changes the chunk, and arda may already have repaired one end of the shipped value."""
    got = anchors.noncanonical(_master([{
        "cdr3.beta": "YLCSSQEGGYGYTF", "__cdr3old.beta": "YLCSSQEGGYGYTFGSG",
        "v.beta": "TRBV29-1*01", "j.beta": "TRBJ1-2*01", "__vcanon.beta": False}]))
    assert got.height == 1
    assert got["cdr3.original"][0] == "YLCSSQEGGYGYTFGSG"
    assert got["repair"][0] == "CSSQEGGYGYTF"


def test_a_chain_with_no_submitted_value_falls_back_to_the_shipped_one():
    """`__cdr3old` is null when nothing needed fixing, in which case the two are the same value."""
    got = anchors.noncanonical(_master([{
        "cdr3.beta": "CASSNEKLF", "__cdr3old.beta": None,
        "v.beta": "TRBV6-1*01", "j.beta": "TRBJ1-4*01"}]))
    assert got.height == 1
    assert got["cdr3.original"][0] == "CASSNEKLF"
    assert got["repair"][0] == "CASSNEKLFF"


def test_nothing_flagged_gives_an_empty_frame_and_an_empty_report():
    """The report is skipped rather than printing a heading with nothing under it."""
    clean = _master([{"cdr3.beta": "CASSNEKLFF", "v.beta": "TRBV6-1*01", "j.beta": "TRBJ1-4*01"}])
    got = anchors.noncanonical(clean)
    assert got.is_empty()
    assert anchors.report(got) == ""


def test_the_report_names_the_defect_and_says_it_blocks_nothing():
    got = anchors.noncanonical(_master([
        {"cdr3.beta": "CASSNEKLF", "v.beta": "TRBV6-1*01", "j.beta": "TRBJ1-4*01"},
        {"cdr3.alpha": "CPDYANKMIF", "v.alpha": "TRAV8D-1*01", "j.alpha": "TRAJ47*01",
         "species": "MusMusculus"}]))
    text = anchors.report(got)
    assert "blocks nothing" in text
    assert "J absent anchor" in text
    assert "anchor table suspect" in text
    assert "CASSNEKLFF" in text


def test_an_empty_corpus_does_not_raise():
    assert anchors.noncanonical(_master([])).is_empty()


def test_the_anchor_table_being_empty_raises_rather_than_reporting_nothing():
    """`load_anchors` returns an empty dict when arda has no reference, so every junction would come
    back unchecked and the report would look like it ran. Measured: all 990 flagged chains, before
    `ensure_reference` was called here."""
    anchors._anchors.cache_clear()
    import arda.cdr3fix

    original = arda.cdr3fix.load_anchors

    def empty(organism):
        return {}

    # `ensure_reference` calls `load_anchors.cache_clear()`, so the stand-in needs one - without it
    # the test fails on an AttributeError instead of on the behaviour it is about.
    empty.cache_clear = lambda: None
    arda.cdr3fix.load_anchors = empty
    try:
        with pytest.raises(RuntimeError, match=r"anchor table|germline reference"):
            anchors.templated("HomoSapiens", "V", "TRBV6-1*01")
    finally:
        arda.cdr3fix.load_anchors = original
        anchors._anchors.cache_clear()
