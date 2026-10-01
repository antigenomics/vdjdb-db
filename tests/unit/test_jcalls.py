"""J calls the junction contradicts (#681): the rule, and what it refuses."""
from __future__ import annotations

import polars as pl

from vdjdb.curate.jcalls import MIN_RUN, contradicted, report


def test_a_junction_that_ends_in_another_gene_names_it():
    """`CAVSDGNAGGTSYGKLTF` is called TRAJ47 and ends in eleven residues of TRAJ52."""
    got = contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTF", "TRAJ47*01")
    assert got is not None
    best, run, called_run = got
    assert best == "TRAJ52"
    assert run >= MIN_RUN and called_run < run


def test_a_call_its_own_junction_supports_is_not_reported():
    """The common case, and the one a false positive here would drown."""
    assert contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTF", "TRAJ52*01") is None


def test_a_species_with_no_germline_reference_is_declined_rather_than_guessed():
    assert contradicted("RattusNorvegicus", "CAVSDGNAGGTSYGKLTF", "TRAJ47*01") is None


def test_an_empty_call_or_junction_is_not_a_finding():
    assert contradicted("HomoSapiens", "", "TRAJ47*01") is None
    assert contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTF", "") is None


def test_the_anchor_cannot_vote_on_its_own_diagnosis():
    """Both sides are compared with the final residue stripped, so substituting it changes nothing.

    A junction whose anchor is corrupt is exactly the population `curate.anchors` reports, and
    letting that residue into this comparison would make one defect produce two findings. The test
    substitutes the anchor rather than appending to it: an appended residue is a longer sequence, a
    different defect, and one this check is right to decline.
    """
    intact = contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTF", "TRAJ47*01")
    corrupt = contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTX", "TRAJ47*01")
    assert intact is not None and corrupt is not None
    assert intact == corrupt, "the gene the junction names must not depend on its last residue"
    assert contradicted("HomoSapiens", "CAVSDGNAGGTSYGKLTFX", "TRAJ47*01") is None


def test_the_report_joins_back_to_the_chain_and_its_chunk():
    chains = pl.DataFrame({
        "record_id": ["r1", "r2"], "gene": ["TRA", "TRA"],
        "cdr3": ["CAVSDGNAGGTSYGKLTF", "CAVSDGNAGGTSYGKLTF"],
        "j.segm": ["TRAJ47*01", "TRAJ52*01"],
    })
    records = pl.DataFrame({
        "record_id": ["r1", "r2"], "species": ["HomoSapiens", "HomoSapiens"],
        "chunk.file": ["a.txt", "b.txt"],
    })
    got = report(chains, records)
    assert got["record_id"].to_list() == ["r1"], "only the contradicted call is a row"
    assert got["chunk.file"].item() == "a.txt"
    assert got["best.gene"].item() == "TRAJ52"


# --- the shipped rule (`recall`) and where it puts the gene (`place`) ---------------------------

def test_a_j_that_misses_the_last_three_residues_is_replaced_by_the_one_gene_that_fits():
    from vdjdb.curate.jcalls import recall

    # TRAJ22 matches 2 residues of this junction, TRAJ45 matches 3, and nothing else reaches 3.
    assert recall("HomoSapiens", "TRA", "CALGGLTF", "TRAJ22*01") == "TRAJ45"


def test_one_mismatch_beside_the_anchor_does_not_cost_the_call():
    """`ADGLPF` against TRAJ45's `ADGLTF`: the residue beside the anchor is skipped, the rest agrees."""
    from vdjdb.curate.jcalls import recall

    assert recall("HomoSapiens", "TRA", "CAASGGGADGLPF", "TRAJ45*01") is None


def test_a_call_that_matches_three_residues_stands_even_where_another_gene_matches_more():
    """`CAMSGGKLIF` ends `GKLI`: TRAJ37 matches 4, TRAJ23 matches 5. Gaining a residue is not enough."""
    from vdjdb.curate.jcalls import recall

    assert recall("HomoSapiens", "TRA", "CAMSGGKLIF", "TRAJ37*01") is None


def test_ties_are_not_broken(monkeypatch):
    from vdjdb.curate import jcalls

    genes = {"TRAJ1": "QQQQAB", "TRAJ2": "QQQQAB", "TRAJ3": "WWWWWW"}
    monkeypatch.setattr(jcalls, "_bodies", lambda species, locus: genes)
    assert jcalls.recall("HomoSapiens", "TRA", "CXXXQQQQABF", "TRAJ3*01") is None


def test_a_call_nothing_else_explains_stands(monkeypatch):
    from vdjdb.curate import jcalls

    monkeypatch.setattr(jcalls, "_bodies", lambda species, locus: {"TRAJ1": "AAAAAA", "TRAJ2": "BBBBBB"})
    assert jcalls.recall("HomoSapiens", "TRA", "CXXXXXXF", "TRAJ1*01") is None


def test_an_unknown_gene_or_species_is_not_judged():
    from vdjdb.curate.jcalls import recall

    assert recall("HomoSapiens", "TRA", "CALGGLTF", "TRAJ999*01") is None
    assert recall("RattusNorvegicus", "TRA", "CALGGLTF", "TRAJ22*01") is None
    assert recall("HomoSapiens", "TRA", "", "TRAJ22*01") is None


def test_place_names_the_allele_the_start_and_whether_the_junction_closes_on_its_anchor():
    from vdjdb.curate.jcalls import place

    allele, start, canonical = place("HomoSapiens", "CAVPMYSGGGADGLAF", "TRAJ45")
    assert allele.startswith("TRAJ45*") and canonical is True
    assert 0 <= start < len("CAVPMYSGGGADGLAF") - 3


def test_the_markup_ships_the_rules_gene_and_keeps_the_engines_own_call_beside_it():
    from vdjdb.annotate.cdr3fix import markup

    keys = pl.DataFrame({"species": ["HomoSapiens"], "cdr3": ["CAVPMYSGGGADGLAF"],
                         "v": ["TRAV12-2"], "j": ["TRAJ45"]})
    row = markup(keys, "alpha").row(0, named=True)
    assert row["__j"].startswith("TRAJ45*"), "the shipped call is the rule's gene"
    assert not row["__jarda"].startswith("TRAJ45"), "arda's own call stays in j.segm.arda"
    assert row["__jcanon"] is True
    off = markup(keys, "alpha", recall_j=False).row(0, named=True)
    assert off["__j"] == row["__jarda"], "without the rule the engine's call ships"
