"""The anchor repair must be no worse than the k-mer scanner it replaced.

`arda.cdr3fix` leads the retired scanner on coordinates given a correct call - measured against
nucleotide truth, `j.start` exact on 8,164 records against 7,987. It does not lead on *tolerating* a
junction that carries framework past an anchor: it trims one end and leaves the other, and
`vdjdb.curate.anchors` computed the repair for the other end and only reported it. The retired build
applied it, so the shipped junction regressed.

Measured on the 2026-09-30 corpus, applying the proposals ungated was worse than doing nothing: 191 of
218 distinct repairs produced a sequence the 2026-06-03 release never ships, because `classify` takes
the first position where the alignment improves and on a junction whose last residue is *mis-read*
that is an internal anchor residue. Gated on germline depth at the J end, 33 repairs apply: 17 are the
release's own value, 16 are more correct than the release, and none is wrong.

The 16 are the per-segment anchor cases. `TRAJ35*01` templates `IGFGNVLHC` and mouse `TRAJ7*01`
templates `DYSNNRLTL`, so a correct junction on either ends in Cys or Leu - and the release's fixed
"ends in Phe or Trp" test wrote an F there. Agreement against the release therefore *falls* by 16 on
these, and that is the repair working.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate.anchors import repairs, templated

#: ``(submitted, species, v, j, expected repair or None)``. The first four are the framework-flanked
#: junctions the release repairs and this build did not: `PMID_11046006` submits
#: ``YFCASSQSPGGVAFFGQG`` and the release ships ``CASSQSPGGVAFF``.
CASES: tuple[tuple[str, str, str, str, str | None], ...] = (
    ("YFCASSQSPGGVAFFGQG", "HomoSapiens", "TRBV14*01", "TRBJ1-1*01", "CASSQSPGGVAFF"),
    ("YFCAVVGTGLGYTFGSG", "HomoSapiens", "TRBV9*01", "TRBJ1-2*01", "CAVVGTGLGYTF"),
    ("YLCSSQEGGYGYTFGSG", "HomoSapiens", "TRBV29-1*01", "TRBJ1-2*01", "CSSQEGGYGYTF"),
    ("YICSDDKYGNTIYFGEG", "HomoSapiens", "TRBV20-1*01", "TRBJ1-3*01", "CSDDKYGNTIYF"),
    # Refused: the terminal residue is mis-read, not followed by framework, so the only anchor
    # residue available is an internal one. Trimming to it destroys the junction - the release ships
    # `CAISGEFGSGANVLTF`, and the ungated pass produced `CAISGEF`.
    ("CAISGEFGSGA", "HomoSapiens", "TRBV10-1*01", "TRBJ2-6*01", None),
    ("CSAGRDGTNEKLFL", "HomoSapiens", "TRBV29-1*01", "TRBJ1-4*01", None),
    ("CASSLGDRAFRNIQY", "HomoSapiens", "TRBV7-9*01", "TRBJ2-4*01", None),
    ("CASSDNPLVGGFTDTQY", "HomoSapiens", "TRBV9*01", "TRBJ2-3*01", None),
    ("CATILSTGGFKTIS", "HomoSapiens", "TRAV17*01", "TRAJ9*01", None),
)


def _repair(submitted: str, species: str, v: str, j: str) -> str | None:
    keys = pl.DataFrame({"species": [species], "cdr3": [submitted], "v": [v], "j": [j]})
    got = repairs(keys)
    return None if got.is_empty() else got["__repaired"].item()


@pytest.mark.parametrize(("submitted", "species", "v", "j", "want"), CASES)
def test_the_repair_matches_the_release_or_declines(submitted, species, v, j, want) -> None:
    assert _repair(submitted, species, v, j) == want


def test_a_repair_never_ends_on_a_residue_the_j_germline_does_not_name() -> None:
    """The gate, stated directly: this is what the ungated pass violated 191 times."""
    for submitted, species, v, j, _ in CASES:
        got = _repair(submitted, species, v, j)
        if got is None:
            continue
        germline = templated(species, "J", j)
        assert germline, f"{j} has no germline, so no repair should have been proposed"
        assert got.endswith(germline[-1]), f"{submitted} -> {got} does not end on {j}'s anchor"
        assert got.startswith("C"), f"{submitted} -> {got} does not start at Cys104"


def test_the_anchor_is_read_per_segment_and_is_not_always_phe_or_trp() -> None:
    """`TRAJ35*01` ends in Cys and mouse `TRAJ7*01` in Leu.

    These are the 16 repairs that disagree with the 2026-06-03 release *because the release is wrong*:
    its fixed Phe-or-Trp test wrote an F where the germline encodes C or L. A test asserting
    "ends in F or W" would forbid the correct answer, which is why there is no such test here.
    """
    assert templated("HomoSapiens", "J", "TRAJ35*01") == "IGFGNVLHC"
    assert templated("MusMusculus", "J", "TRAJ7*01") == "DYSNNRLTL"
    assert _repair("CALEGFGNVLHF", "HomoSapiens", "TRAV21*02", "TRAJ35*01") == "CALEGFGNVLHC"
    assert _repair("CAVSIDYSNNRLTLF", "MusMusculus", "TRAV9D-4*01", "TRAJ7*01") == "CAVSIDYSNNRLTL"


def test_a_call_fix_is_not_turned_into_a_sequence_fix() -> None:
    """Where a functional sibling allele explains the junction, the call is wrong and the sequence is
    right. Rewriting the residue would destroy the evidence for the real defect - the reasoning
    `MAX_REPLACE = 0` already applies in `vdjdb.annotate.cdr3fix`. Mouse `TRAJ47` resolves to `*01`,
    an ORF allele templating `HYANKMIC`, while every corpus chain reads `DYANKMIF`, which is `*02`.
    """
    assert _repair("CAAGDYANKMIF", "MusMusculus", "TRAV16D/DV11*01", "TRAJ47*01") is None
