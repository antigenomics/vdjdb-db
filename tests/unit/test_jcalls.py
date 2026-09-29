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
