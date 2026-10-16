"""Synthetic checks for the dCODE chain graph. No cutoff grid."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
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


def _singleton_panel(
    n_background: int,
    extras: list[tuple[str, dict[str, float], float]],
    columns: list[str],
    nc_col: str,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Background labels at 0, plus extra labels of (clone_id, dextramer UMIs, max NC)."""
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(n_background):
        count_rows.append({"barcode": f"z{i}", **{col: 0.0 for col in columns}})
        control_rows.append({"barcode": f"z{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"z{i}", "clone_id": f"z{i}"})
    for clone_id, dex, nc in extras:
        barcode = f"{clone_id}-0"
        count_rows.append({"barcode": barcode, **{col: float(dex.get(col, 0.0)) for col in columns}})
        control_rows.append({"barcode": barcode, nc_col: float(nc)})
        clone_rows.append({"donor": 1, "barcode": barcode, "clone_id": clone_id})
    return pl.DataFrame(clone_rows), pl.DataFrame(count_rows), pl.DataFrame(control_rows)


def test_control_edge_is_mean_plus_three_sd_not_the_maximum():
    """One control cell at 50 UMIs must not set the bar at 50."""
    nc_max = np.array([0.0] * 42 + [50.0])
    mean, sd, fence = dcode.control_edge(nc_max)
    assert mean < 0.2
    assert fence < math.log1p(11)
    assert fence > math.log1p(1)
    assert fence < math.log1p(50)
    assert sd > 0


def test_ela_singleton_clears_its_own_control_and_not_the_donor_maximum():
    """ELAGIGILTV at 11 dextramer UMIs, own NC 0, is a call. A count of 1 is not.

    The donor's highest negative-control cell is 50. The score is log1p of the
    dextramer UMI, and n_cells of the hit is 1.
    """
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        40,
        [
            ("donor_max", {ela: 0.0}, 50.0),
            ("noise", {ela: 1.0}, 0.0),
            ("ela", {ela: 11.0}, 0.0),
        ],
        [ela],
        nc_col,
    )
    # The helper stores the NC on the control frame via the third field, and
    # also copies a per-cell override when present. donor_max's NC is 50.
    calls = dcode.call_specificities(clones, counts, [ela], controls, [nc_col])
    assert calls["clone_id"].to_list() == ["ela"]
    assert calls["peptide"].to_list() == ["ELAGIGILTV"]
    assert calls["mhc"].to_list() == ["HLA-A*02:01"]
    assert calls["n_cells"].to_list() == [1]
    assert math.isclose(float(calls["score"][0]), math.log1p(11), rel_tol=0, abs_tol=1e-9)


def test_sticky_peptide_high_on_every_label_is_not_a_call():
    """KLGGALQAK at 20 on every label beats a zero control and is still not a call.

    ELAGIGILTV at 40 on one label, against a bulk of zeros, is a call.
    """
    klg = "A0301_KLGGALQAK_IE-1_CMV"
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(30):
        count_rows.append({"barcode": f"b{i}", klg: 20.0, ela: 0.0})
        control_rows.append({"barcode": f"b{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"b{i}", "clone_id": f"c{i}"})
    for i in range(4):
        count_rows.append({"barcode": f"m{i}", klg: 20.0, ela: 40.0})
        control_rows.append({"barcode": f"m{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"m{i}", "clone_id": "mart"})
    counts = pl.DataFrame(count_rows)
    controls = pl.DataFrame(control_rows)
    dex = [klg, ela]
    assert dcode.dextramer_columns([*dex, nc_col]) == dex
    assert dcode.negative_control_columns([*dex, nc_col]) == [nc_col]
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows), counts, dex, controls, [nc_col],
    )
    assert set(calls["peptide"].to_list()) == {"ELAGIGILTV"}
    assert calls["clone_id"].to_list() == ["mart"]
    assert calls["mhc"].to_list() == ["HLA-A*02:01"]
    assert int(calls["n_cells"][0]) == 4
    assert math.isclose(float(calls["score"][0]), math.log1p(40), rel_tol=0, abs_tol=1e-9)


def test_eight_cells_at_one_umi_are_not_one_cell_at_eight():
    """Cell count is the sample size. Eight dextramer UMIs of 1 do not sum to 8."""
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(40):
        count_rows.append({"barcode": f"z{i}", pep: 0.0})
        control_rows.append({"barcode": f"z{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"z{i}", "clone_id": f"z{i}"})
    for i in range(4):
        count_rows.append({"barcode": f"n{i}", pep: 0.0})
        control_rows.append({"barcode": f"n{i}", nc_col: 3.0})
        clone_rows.append({"donor": 1, "barcode": f"n{i}", "clone_id": f"n{i}"})
    for i in range(8):
        count_rows.append({"barcode": f"e{i}", pep: 1.0})
        control_rows.append({"barcode": f"e{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"e{i}", "clone_id": "eight"})
    count_rows.append({"barcode": "one", pep: 8.0})
    control_rows.append({"barcode": "one", nc_col: 0.0})
    clone_rows.append({"donor": 1, "barcode": "one", "clone_id": "one"})
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows),
        pl.DataFrame(count_rows),
        [pep],
        pl.DataFrame(control_rows),
        [nc_col],
    )
    assert calls["clone_id"].to_list() == ["one"]
    assert calls["n_cells"].to_list() == [1]
    assert math.isclose(float(calls["score"][0]), math.log1p(8), rel_tol=0, abs_tol=1e-9)
    assert float(calls["score"][0]) != 8.0


def test_within_label_median_ignores_one_hot_cell():
    """One cell at 100 dextramer UMIs among three zeros has median log1p 0.

    Four cells stable at 30 are the within-label signal. The raw mean of the
    hot cell would have been 25 and is not the score.
    """
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(20):
        count_rows.append({"barcode": f"z{i}", pep: 0.0})
        control_rows.append({"barcode": f"z{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"z{i}", "clone_id": f"z{i}"})
    for i, value in enumerate([0.0, 0.0, 0.0, 100.0]):
        count_rows.append({"barcode": f"h{i}", pep: value})
        control_rows.append({"barcode": f"h{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"h{i}", "clone_id": "hot"})
    for i in range(4):
        count_rows.append({"barcode": f"s{i}", pep: 30.0})
        control_rows.append({"barcode": f"s{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"s{i}", "clone_id": "stable"})
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows),
        pl.DataFrame(count_rows),
        [pep],
        pl.DataFrame(control_rows),
        [nc_col],
    )
    assert calls["clone_id"].to_list() == ["stable"]
    assert calls["n_cells"].to_list() == [4]
    assert math.isclose(float(calls["score"][0]), math.log1p(30), rel_tol=0, abs_tol=1e-9)


def test_cross_reactive_label_keeps_both_peptides():
    """Two peptides that each clear the control and the across-label bulk both stay."""
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    gil = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nlv = "A0201_NLVPMVATV_pp65_CMV"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        40,
        [
            ("donor_max", {ela: 0.0, gil: 0.0, nlv: 0.0}, 50.0),
            ("both", {ela: 40.0, gil: 40.0, nlv: 1.0}, 0.0),
        ],
        [ela, gil, nlv],
        nc_col,
    )
    calls = dcode.call_specificities(clones, counts, [ela, gil, nlv], controls, [nc_col])
    assert set(calls["peptide"].to_list()) == {"ELAGIGILTV", "GILGFVFTL"}
    assert set(calls["clone_id"].to_list()) == {"both"}


def test_contig_umi_column_is_not_the_dextramer_score():
    """A contig-like UMI column of 500 is not read. The score is log1p(11)."""
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        40,
        [
            ("donor_max", {ela: 0.0}, 50.0),
            ("ela", {ela: 11.0}, 0.0),
        ],
        [ela],
        nc_col,
    )
    counts = counts.with_columns(pl.lit(500.0).alias("umis"))
    calls = dcode.call_specificities(clones, counts, [ela], controls, [nc_col])
    assert calls.height == 1
    assert math.isclose(float(calls["score"][0]), math.log1p(11), rel_tol=0, abs_tol=1e-9)


def test_orphan_alpha_and_beta_are_labels_and_a_second_beta_is_not():
    """A barcode with one chain is a label. Two betas are not."""
    rows = [
        _contig("pair", "TRA", "1", "CAAA", 10),
        _contig("pair", "TRB", "1", "CASS", 10),
        _contig("betaonly", "TRB", "2", "CASSB", 12),
        _contig("alphaonly", "TRA", "3", "CAAAA", 12),
        _contig("doublet", "TRA", "4", "CADD", 10),
        _contig("doublet", "TRB", "5", "CASS1", 12),
        _contig("doublet", "TRB", "6", "CASS2", 11),
    ]
    clones = dcode.clone_table(dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1)))
    shapes = dict(zip(clones["barcode"].to_list(), clones["shape"].to_list(), strict=True))
    assert shapes["pair"] == "ab"
    assert shapes["betaonly"] == "b"
    assert shapes["alphaonly"] == "a"
    assert "doublet" not in shapes
    beta_id = clones.filter(pl.col("barcode") == "betaonly")["clone_id"][0]
    alpha_id = clones.filter(pl.col("barcode") == "alphaonly")["clone_id"][0]
    assert beta_id.startswith("||TRB|")
    assert alpha_id.endswith("||")
    assert "TRB|" not in alpha_id.split("||", 1)[1]

    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    panel_clones, counts, controls = _singleton_panel(
        40,
        [("donor_max", {ela: 0.0}, 50.0)],
        [ela],
        nc_col,
    )
    extra_clones = pl.DataFrame([
        {"donor": 1, "barcode": "betaonly", "clone_id": beta_id},
    ])
    extra_counts = pl.DataFrame([{"barcode": "betaonly", ela: 11.0}])
    extra_controls = pl.DataFrame([{"barcode": "betaonly", nc_col: 0.0}])
    calls = dcode.call_specificities(
        pl.concat([panel_clones, extra_clones], how="diagonal_relaxed"),
        pl.concat([counts, extra_counts], how="diagonal_relaxed"),
        [ela],
        pl.concat([controls, extra_controls], how="diagonal_relaxed"),
        [nc_col],
    )
    assert calls["clone_id"].to_list() == [beta_id]
    assert calls["n_cells"].to_list() == [1]


def test_beta_only_dextramer_umi_is_not_rescaled_by_contig_umi(tmp_path: Path) -> None:
    """A leaked chain's contig UMI must not pull down the dextramer median.

    The beta-only cells carry the binder beta plus a low contig-UMI copy of the
    paired clone's alpha. The score is log1p of the measured dextramer UMI.
    """
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    rows = []
    features = []
    for i in range(30):
        rows.append(_contig(f"z{i}", "TRA", str(i), f"CAV{i}", 12))
        rows.append(_contig(f"z{i}", "TRB", str(i), f"CASS{i}", 12))
        features.append({"barcode": f"z{i}", pep: 0.0, nc_col: 0.0})
    for i in range(5):
        rows.append(_contig(f"p{i}", "TRA", "G", "CAAA", 30))
        rows.append(_contig(f"p{i}", "TRB", "G", "CASSGILG", 40))
        features.append({"barcode": f"p{i}", pep: 100.0, nc_col: 0.0})
    for i in range(4):
        rows.append(_contig(f"b{i}", "TRB", "G", "CASSGILG", 40))
        rows.append(_contig(f"b{i}", "TRA", "G", "CAAA", 1))
        features.append({"barcode": f"b{i}", pep: 100.0, nc_col: 0.0})
    features.append({"barcode": "n0", pep: 0.0, nc_col: 20.0})
    rows.append(_contig("n0", "TRA", "N", "CANEG", 10))
    rows.append(_contig("n0", "TRB", "N", "CASSNEG", 10))
    donor_dir = tmp_path
    pl.DataFrame(rows).write_csv(
        donor_dir / "vdj_v1_hs_aggregated_donor1_all_contig_annotations.csv.gz",
        compression="gzip",
    )
    pl.DataFrame(features).write_csv(
        donor_dir / "vdj_v1_hs_aggregated_donor1_binarized_matrix.csv.gz",
        compression="gzip",
    )
    result = dcode.process_donor(donor_dir, 1)
    beta_only = result["clones"].filter(pl.col("shape") == "b")
    assert beta_only.height == 4
    hit = result["calls"].filter(pl.col("clone_id") == beta_only["clone_id"][0])
    assert hit.height == 1
    assert hit["peptide"].to_list() == ["GILGFVFTL"]
    assert math.isclose(float(hit["score"][0]), math.log1p(100), rel_tol=0, abs_tol=1e-9)
