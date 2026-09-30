"""The ``epitopes`` catalogue: one row per (epitope, species), with what supports it.

The proteome link (#632) is read from ``proofreading/epitope_proteome.tsv``, a committed reviewed
input, never resolved during a build - `mhcmatch` fetches the proteome from HuggingFace, so doing it
here would put a network call in the critical path and make the answer depend on the day.
"""
from __future__ import annotations

import polars as pl

from vdjdb.assemble.epitopes import build_epitopes


def test_the_proteome_link_is_read_from_the_committed_authority(tmp_path):
    """#632. `SLLMWITQV` is 29,729 records and `SLLMWITQC` is 13, and nothing linked them.

    Not a defect flag: the epitope sequence is the ground truth - it is the peptide the experiment
    used - and the difference from the proteome is almost always deliberate. The column exists so a
    query for one form can reach the other.
    """
    (tmp_path / "proofreading").mkdir()
    (tmp_path / "proofreading" / "epitope_proteome.tsv").write_text(
        "antigen.epitope\tantigen.species\tantigen.gene\tverdict\tsource.protein\tsource.gene\t"
        "source.position\tsource.peptide\tsource.subs\trecords\treferences\n"
        "SLLMWITQV\tHomoSapiens\tNY-ESO-1\tone_substitution\tCTG1B_HUMAN\tCTAG1A\t157\t"
        "SLLMWITQC\t9C>V\t29729\t41\n"
        "SLLMWITQC\tHomoSapiens\tNY-ESO-1\texact\tCTG1B_HUMAN\tCTAG1A\t157\t\t\t13\t2\n")
    records = pl.DataFrame({
        "antigen.epitope": ["SLLMWITQV", "SLLMWITQC"],
        "antigen.species": ["HomoSapiens", "HomoSapiens"],
        "antigen.gene": ["NY-ESO-1", "NY-ESO-1"],
        "mhc.class": ["MHCI", "MHCI"], "record_id": ["VDJDB1", "VDJDB2"],
        "reference.id": ["PMID:1", "PMID:2"]})
    chains = pl.DataFrame({"record_id": ["VDJDB1", "VDJDB2"], "clonotype_id": ["CT1", "CT2"]})

    got = build_epitopes(records, chains, root=tmp_path).sort("antigen.epitope")
    assert got["proteome.peptide"].to_list() == ["", "SLLMWITQC"]
    assert got["proteome.substitution"].to_list() == ["", "9C>V"]


def test_an_epitope_the_authority_does_not_cover_gets_an_empty_link(tmp_path):
    """Most rows, and every viral or bacterial epitope: only the two self proteomes are read.

    Empty string, not null - hard rule 6.
    """
    (tmp_path / "proofreading").mkdir()
    (tmp_path / "proofreading" / "epitope_proteome.tsv").write_text(
        "antigen.epitope\tantigen.species\tantigen.gene\tverdict\tsource.protein\tsource.gene\t"
        "source.position\tsource.peptide\tsource.subs\trecords\treferences\n")
    records = pl.DataFrame({
        "antigen.epitope": ["GILGFVFTL"], "antigen.species": ["InfluenzaA"],
        "antigen.gene": ["M"], "mhc.class": ["MHCI"], "record_id": ["VDJDB1"],
        "reference.id": ["PMID:1"]})
    chains = pl.DataFrame({"record_id": ["VDJDB1"], "clonotype_id": ["CT1"]})
    got = build_epitopes(records, chains, root=tmp_path)
    assert got["proteome.peptide"].to_list() == [""]
    assert got["proteome.substitution"].to_list() == [""]


def test_one_peptide_in_two_proteomes_states_the_link_once(tmp_path):
    """`VEALYLVSG` has a row per proteome - human `INS` and mouse `Ins2` both answer `VEALYLVCG`
    at `8C>S` - which is the same statement twice rather than two, so the join must not duplicate
    the epitope row."""
    (tmp_path / "proofreading").mkdir()
    (tmp_path / "proofreading" / "epitope_proteome.tsv").write_text(
        "antigen.epitope\tantigen.species\tantigen.gene\tverdict\tsource.protein\tsource.gene\t"
        "source.position\tsource.peptide\tsource.subs\trecords\treferences\n"
        "VEALYLVSG\tHomoSapiens\tINS\tone_substitution\tINS_HUMAN\tINS\t40\tVEALYLVCG\t8C>S\t2495\t9\n"
        "VEALYLVSG\tMusMusculus\tIns2\tone_substitution\tINS2_MOUSE\tIns1\t40\tVEALYLVCG\t8C>S\t94\t3\n")
    records = pl.DataFrame({
        "antigen.epitope": ["VEALYLVSG"], "antigen.species": ["HomoSapiens"],
        "antigen.gene": ["INS"], "mhc.class": ["MHCI"], "record_id": ["VDJDB1"],
        "reference.id": ["PMID:1"]})
    chains = pl.DataFrame({"record_id": ["VDJDB1"], "clonotype_id": ["CT1"]})
    got = build_epitopes(records, chains, root=tmp_path)
    assert got.height == 1, "the epitope row must not fan out over the authority's rows"
    assert got["proteome.peptide"][0] == "VEALYLVCG"
