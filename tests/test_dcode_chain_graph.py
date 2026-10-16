"""Synthetic checks for the dCODE chain graph. No cutoff grid."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import polars as pl

_SCRIPT = Path(__file__).resolve().parents[1] / "processing" / "dcode_chain_graph.py"
_spec = importlib.util.spec_from_file_location("dcode_chain_graph", _SCRIPT)
assert _spec is not None and _spec.loader is not None
dcode = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dcode)


def _contig(barcode: str, locus: str, gene: str, cdr3: str, umi: int, *, cell: str = "True") -> dict:
    return {
        "barcode": barcode,
        "is_cell": cell,
        "chain": locus,
        "v_gene": f"{locus}V-{gene}",
        "j_gene": f"{locus}J-{gene}",
        "productive": "True",
        "cdr3": cdr3,
        "umis": umi,
    }


def test_boolean_is_cell_from_polars_csv_inference():
    contigs = pl.DataFrame({
        "barcode": ["c1"],
        "is_cell": [True],
        "chain": ["TRA"],
        "v_gene": ["TRAV1"],
        "j_gene": ["TRAJ1"],
        "productive": ["True"],
        "cdr3": ["CAAA"],
        "umis": [4],
    })
    chains = dcode.chain_table(contigs, donor=1)
    assert chains["is_cell"][0] is True


def test_boolean_productive_is_kept():
    """Polars infers ``is_cell`` as Boolean. A Boolean ``productive`` must not empty the table."""
    contigs = pl.DataFrame({
        "barcode": ["c1", "c1", "c2"],
        "is_cell": [True, True, False],
        "chain": ["TRA", "TRB", "TRA"],
        "v_gene": ["TRAV1", "TRBV1", "TRAV2"],
        "j_gene": ["TRAJ1", "TRBJ1", "TRAJ2"],
        "productive": [True, False, True],
        "cdr3": ["CAAA", "CASS", "CAGG"],
        "umis": [4, 8, 3],
    })
    chains = dcode.chain_table(contigs, donor=1)
    assert chains.height == 2
    assert set(chains.filter(pl.col("is_cell"))["cdr3"].to_list()) == {"CAAA"}


def test_lowercase_true_is_a_cell_and_edges_are_undirected():
    contigs = pl.DataFrame([
        _contig("c1", "TRA", "1", "CAAA", 5),
        _contig("c1", "TRB", "1", "CASS", 8),
        _contig("c2", "TRA", "2", "CAGG", 4, cell="False"),
    ])
    chains = dcode.chain_table(contigs, donor=1)
    assert chains.filter(pl.col("is_cell")).height == 2
    edges = dcode.edge_table(chains)
    assert edges.height == 1
    assert edges["chain_a"][0] < edges["chain_b"][0]


def test_leakage_of_an_expanded_beta_is_not_a_second_clone():
    rows = [
        _contig("true", "TRA", "1", "CAAA", 20),
        _contig("true", "TRB", "9", "CASSLEAK", 50),
    ]
    for i in range(40):
        rows.append(_contig(f"l{i}", "TRA", str(i), f"CAV{i}", 15))
        rows.append(_contig(f"l{i}", "TRB", str(i), f"CASS{i}", 15))
        rows.append(_contig(f"l{i}", "TRB", "9", "CASSLEAK", 1))
    chains = dcode.chain_table(pl.DataFrame(rows), donor=1)
    endogenous = dcode.mark_endogenous(chains)
    leak_cells = endogenous.filter(pl.col("barcode") != "true")
    leaked = leak_cells.filter(pl.col("cdr3") == "CASSLEAK")
    assert leaked.height == 40
    assert int(leaked["endogenous"].sum()) == 0
    own = leak_cells.filter(pl.col("cdr3") != "CASSLEAK")
    assert int(own["endogenous"].sum()) == 80
    true_beta = endogenous.filter((pl.col("barcode") == "true") & (pl.col("cdr3") == "CASSLEAK"))
    assert bool(true_beta["endogenous"][0])
    clones = dcode.clone_table(endogenous)
    assert clones.filter(pl.col("barcode") == "true").height == 1
    assert clones.filter(pl.col("barcode") == "true")["shape"][0] == "ab"
    assert clones.height == 41


def test_repeated_dual_alpha_is_kept():
    rows = []
    for i in range(5):
        rows.append(_contig(f"d{i}", "TRA", "1", "CAAA", 20))
        rows.append(_contig(f"d{i}", "TRA", "2", "CABB", 15))
        rows.append(_contig(f"d{i}", "TRB", "1", "CASS", 20))
    endogenous = dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1))
    clones = dcode.clone_table(endogenous)
    assert clones.height == 5
    assert set(clones["shape"].to_list()) == {"aab"}


def test_second_beta_is_not_a_clone():
    rows = [
        _contig("x", "TRA", "1", "CAAA", 10),
        _contig("x", "TRB", "1", "CASS1", 12),
        _contig("x", "TRB", "2", "CASS2", 11),
    ]
    endogenous = dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1))
    # Both betas clear a zero ambient rate, so the barcode is two betas and is dropped.
    assert dcode.clone_table(endogenous).height == 0


def test_additive_doublet_is_dropped_and_dual_alpha_feature_is_kept():
    def chain(locus: str, gene: str, cdr3: str) -> str:
        return f"{locus}|{locus}V-{gene}|{locus}J-{gene}|{cdr3}"

    a1 = chain("TRA", "1", "CAAA")
    a2 = chain("TRA", "2", "CABB")
    b1 = chain("TRB", "1", "CASS1")
    b2 = chain("TRB", "2", "CASS2")
    clone_rows = []
    feature_rows = []
    for i in range(3):
        clone_rows.append({
            "donor": 1, "barcode": f"p{i}", "shape": "ab",
            "alpha1": a1, "alpha2": None, "beta": b1, "clone_id": f"{a1}||{b1}",
        })
        feature_rows.append({"barcode": f"p{i}", "CD8a": 10.0, "MART": 10.0, "KLG": 0.0, "NC": 1.0})
    for i in range(3):
        clone_rows.append({
            "donor": 1, "barcode": f"q{i}", "shape": "ab",
            "alpha1": a2, "alpha2": None, "beta": b2, "clone_id": f"{a2}||{b2}",
        })
        feature_rows.append({"barcode": f"q{i}", "CD8a": 10.0, "MART": 0.0, "KLG": 10.0, "NC": 1.0})
    # Sum of the two parents.
    clone_rows.append({
        "donor": 1, "barcode": "dbl", "shape": "aab",
        "alpha1": a1, "alpha2": a2, "beta": b1, "clone_id": f"{a1}|{a2}||{b1}",
    })
    feature_rows.append({"barcode": "dbl", "CD8a": 20.0, "MART": 10.0, "KLG": 10.0, "NC": 2.0})
    # Same chains, features of one parent only.
    clone_rows.append({
        "donor": 1, "barcode": "dual", "shape": "aab",
        "alpha1": a1, "alpha2": a2, "beta": b1, "clone_id": f"{a1}|{a2}||{b1}",
    })
    feature_rows.append({"barcode": "dual", "CD8a": 10.0, "MART": 10.0, "KLG": 0.0, "NC": 1.0})
    clones = pl.DataFrame(clone_rows)
    features = pl.DataFrame(feature_rows)
    kept = dcode.drop_doublets(clones, features, ["CD8a", "MART", "KLG", "NC"])
    assert "dbl" not in kept["barcode"].to_list()
    assert "dual" in kept["barcode"].to_list()


def test_sticky_column_is_the_negative_control_median_and_a_rare_peptide_is_called():
    rows = []
    klg = "A0301_KLGGALQAK_IE-1_CMV"
    mart = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc = "A0201_ALIAPVHAV_NC"
    for i in range(30):
        rows.append({"barcode": f"b{i}", klg: 20.0, mart: 0.0, nc: 10.0})
    for i in range(4):
        rows.append({"barcode": f"m{i}", klg: 20.0, mart: 40.0, nc: 10.0})
    features = pl.DataFrame(rows)
    dex = dcode.dextramer_columns(features.columns)
    nc = dcode.negative_control_columns(features.columns)
    assert dex == ["A0301_KLGGALQAK_IE-1_CMV", "A0201_ELAGIGILTV_MART-1_Cancer"]
    residuals, gammas = dcode.column_residuals(features, dex, nc)
    assert gammas["A0301_KLGGALQAK_IE-1_CMV"] == 2.0
    assert gammas["A0201_ELAGIGILTV_MART-1_Cancer"] == 0.0
    assert float(residuals.filter(pl.col("barcode") == "m0")["A0301_KLGGALQAK_IE-1_CMV"][0]) == 0.0
    assert float(residuals.filter(pl.col("barcode") == "m0")["A0201_ELAGIGILTV_MART-1_Cancer"][0]) == 40.0
    clones = pl.DataFrame(
        [{"donor": 1, "barcode": f"b{i}", "clone_id": f"c{i}"} for i in range(30)]
        + [{"donor": 1, "barcode": f"m{i}", "clone_id": "mart"} for i in range(4)]
    )
    calls = dcode.call_specificities(clones, residuals, dex, features, nc)
    peptides = set(calls["peptide"].to_list())
    assert "ELAGIGILTV" in peptides
    assert "KLGGALQAK" not in peptides
    assert calls.filter(pl.col("peptide") == "ELAGIGILTV")["mhc"][0] == "HLA-A*02:01"


def test_unit_residual_does_not_clear_the_negative_control_envelope():
    """A pile of residual 1 is inside the control sample. A separated residual is not."""
    mart = "A0201_ELAGIGILTV_MART-1_Cancer"
    vteh = "A0101_VTEHDTLLY_pp65_CMV"
    nc_col = "A0201_ALIAPVHAV_NC"
    residual_rows = []
    control_rows = []
    clone_rows = []
    for i in range(50):
        residual_rows.append({"barcode": f"b{i}", mart: float(i % 2), vteh: 1.0})
        control_rows.append({"barcode": f"b{i}", nc_col: 0.0 if i < 40 else 30.0})
        clone_rows.append({"donor": 1, "barcode": f"b{i}", "clone_id": f"c{i}"})
    residual_rows.append({"barcode": "hi", mart: 80.0, vteh: 1.0})
    control_rows.append({"barcode": "hi", nc_col: 0.0})
    clone_rows.append({"donor": 1, "barcode": "hi", "clone_id": "hi"})
    residuals = pl.DataFrame(residual_rows)
    controls = pl.DataFrame(control_rows)
    clones = pl.DataFrame(clone_rows)
    dex = [mart, vteh]
    calls = dcode.call_specificities(clones, residuals, dex, controls, [nc_col])
    peptides = set(calls["peptide"].to_list()) if calls.height else set()
    assert peptides == {"ELAGIGILTV"}
    assert calls["clone_id"].to_list() == ["hi"]


def test_clone_size_does_not_sum_a_shoulder_into_a_call():
    """Eight cells at residual 6 stay under a control cell of 20. A singleton at 80 clears it."""
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    residual_rows = [{"barcode": f"e{i}", pep: 6.0} for i in range(8)]
    control_rows = [{"barcode": f"e{i}", nc_col: 0.0} for i in range(8)]
    clone_rows = [{"donor": 1, "barcode": f"e{i}", "clone_id": "expanded"} for i in range(8)]
    residual_rows.append({"barcode": "ctrl", pep: 0.0})
    control_rows.append({"barcode": "ctrl", nc_col: 20.0})
    clone_rows.append({"donor": 1, "barcode": "ctrl", "clone_id": "ctrl"})
    residual_rows.append({"barcode": "real", pep: 80.0})
    control_rows.append({"barcode": "real", nc_col: 0.0})
    clone_rows.append({"donor": 1, "barcode": "real", "clone_id": "real"})
    residual_rows.append({"barcode": "tied", pep: 20.0})
    control_rows.append({"barcode": "tied", nc_col: 0.0})
    clone_rows.append({"donor": 1, "barcode": "tied", "clone_id": "tied"})
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows),
        pl.DataFrame(residual_rows),
        [pep],
        pl.DataFrame(control_rows),
        [nc_col],
    )
    assert calls["clone_id"].to_list() == ["real"]


def test_a_mode_below_the_controls_does_not_leak_past_them():
    """A declining tail past the controls is not a call. A mode above them is."""
    sticky = "A0301_KLGGALQAK_IE-1_CMV"
    flu = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    residual_rows = []
    control_rows = []
    clone_rows = []

    def add(barcode: str, klg: float, gilg: float, nc: float) -> None:
        residual_rows.append({"barcode": barcode, sticky: klg, flu: gilg})
        control_rows.append({"barcode": barcode, nc_col: nc})
        clone_rows.append({"donor": 1, "barcode": barcode, "clone_id": barcode})

    for i, value in enumerate([4.0] * 20 + [8.0] * 8 + [16.0] * 4 + [40.0] * 2):
        add(f"s{i}", value, 0.0, 0.0)
    for i, value in enumerate([0.0] * 15 + [4.0] * 2 + [40.0] * 8 + [80.0] * 5):
        add(f"g{i}", 0.0, value, 0.0)
    add("ctrl", 0.0, 0.0, 30.0)
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows),
        pl.DataFrame(residual_rows),
        [sticky, flu],
        pl.DataFrame(control_rows),
        [nc_col],
    )
    assert set(calls["peptide"].to_list()) == {"GILGFVFTL"}
    assert calls.height == 13
