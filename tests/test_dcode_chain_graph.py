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


def test_expanded_clone_is_not_the_ambient_profile():
    """A frequent beta expressed above the extra-chain rate stays endogenous.

    Rare second betas would make a frequency-scaled ambient mean larger than
    this clone's own count. The clone is the receptor. Copies at the extra-chain
    rate, behind another beta, are contamination.
    """
    rows = []
    for i in range(40):
        rows.append(_contig(f"e{i}", "TRA", str(i), f"CAE{i}", 5))
        rows.append(_contig(f"e{i}", "TRB", str(i), f"CASE{i}", 6))
        rows.append(_contig(f"e{i}", "TRB", f"x{i}", f"CASSX{i}", 1))
    for i in range(30):
        rows.append(_contig(f"c{i}", "TRA", f"a{i}", f"CAC{i}", 4))
        rows.append(_contig(f"c{i}", "TRB", "9", "CASSPUBLIC", 6))
    for i in range(10):
        rows.append(_contig(f"l{i}", "TRA", f"b{i}", f"CAL{i}", 8))
        rows.append(_contig(f"l{i}", "TRB", f"b{i}", f"CASSL{i}", 10))
        rows.append(_contig(f"l{i}", "TRB", "9", "CASSPUBLIC", 1))
    endogenous = dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1))
    public = endogenous.filter(pl.col("cdr3") == "CASSPUBLIC")
    high = public.filter(pl.col("umi") == 6)
    low = public.filter(pl.col("umi") == 1)
    assert high.height == 30
    assert int(high["endogenous"].sum()) == 30
    assert low.height == 10
    assert int(low["endogenous"].sum()) == 0
    own = endogenous.filter(
        pl.col("barcode").str.starts_with("c") & (pl.col("cdr3") != "CASSPUBLIC")
    )
    assert own.height == 30
    assert int(own["endogenous"].sum()) == 30


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


def _ab_row(barcode: str, alpha: str, beta: str) -> dict:
    return {
        "donor": 1, "barcode": barcode, "shape": "ab",
        "alpha1": alpha, "alpha2": None, "beta": beta, "clone_id": f"{alpha}||{beta}",
    }


def test_two_clones_of_one_alpha_are_not_a_doublet():
    """Other betas of a public alpha are not two parents of a dual-alpha cell.

    The second alpha has no clone of its own. The peptide count is one spike.
    The barcode stays, because no cross-alpha parent pair exists.
    """
    def chain(locus: str, gene: str, cdr3: str) -> str:
        return f"{locus}|{locus}V-{gene}|{locus}J-{gene}|{cdr3}"

    public = chain("TRA", "12-2", "CAVNN")
    second = chain("TRA", "2", "CAVVV")
    beta = chain("TRB", "4", "CASRM")
    left = chain("TRB", "10", "CAISE")
    right = chain("TRB", "2", "CASSE")
    rows = []
    features = []
    for i in range(4):
        rows.append(_ab_row(f"L{i}", public, left))
        features.append({"barcode": f"L{i}", "ELA": 16.0, "NC": 0.0})
        rows.append(_ab_row(f"R{i}", public, right))
        features.append({"barcode": f"R{i}", "ELA": 7.0, "NC": 0.0})
    rows.append({
        "donor": 1, "barcode": "cell", "shape": "aab",
        "alpha1": public, "alpha2": second, "beta": beta,
        "clone_id": f"{public}|{second}||{beta}",
    })
    features.append({"barcode": "cell", "ELA": 65.0, "NC": 0.0})
    kept = dcode.drop_doublets(pl.DataFrame(rows), pl.DataFrame(features), ["ELA", "NC"])
    assert "cell" in kept["barcode"].to_list()


def test_public_alphas_drop_only_when_counts_fit_the_sum():
    """Two public alphas drop when the sum of their clones fits, and stay when it does not."""
    def chain(locus: str, gene: str, cdr3: str) -> str:
        return f"{locus}|{locus}V-{gene}|{locus}J-{gene}|{cdr3}"

    a1 = chain("TRA", "1", "CAAA")
    a2 = chain("TRA", "2", "CABB")
    beta = chain("TRB", "9", "CASSNEW")
    b1 = chain("TRB", "1", "CASS1")
    b2 = chain("TRB", "2", "CASS2")
    rows = []
    features = []
    for i in range(3):
        rows.append(_ab_row(f"p{i}", a1, b1))
        features.append({"barcode": f"p{i}", "ELA": 10.0, "GIL": 0.0})
        rows.append(_ab_row(f"q{i}", a2, b2))
        features.append({"barcode": f"q{i}", "ELA": 0.0, "GIL": 10.0})
    rows.append({
        "donor": 1, "barcode": "spike", "shape": "aab",
        "alpha1": a1, "alpha2": a2, "beta": beta,
        "clone_id": f"{a1}|{a2}||{beta}",
    })
    features.append({"barcode": "spike", "ELA": 40.0, "GIL": 0.0})
    rows.append({
        "donor": 1, "barcode": "sum", "shape": "aab",
        "alpha1": a1, "alpha2": a2, "beta": beta,
        "clone_id": f"{a1}|{a2}||{beta}",
    })
    features.append({"barcode": "sum", "ELA": 10.0, "GIL": 10.0})
    kept = dcode.drop_doublets(pl.DataFrame(rows), pl.DataFrame(features), ["ELA", "GIL"])
    barcodes = kept["barcode"].to_list()
    assert "spike" in barcodes
    assert "sum" not in barcodes


def test_primary_low_count_stays_endogenous_and_a_later_copy_does_not():
    """A first-slot count of 1 from a receptor chain stays. The same count behind another chain does not.

    The slot prior is estimated from these rows. It is not a count cutoff.
    """
    rows = []
    for i in range(30):
        rows.append(_contig(f"e{i}", "TRA", str(i), f"CAE{i}", 6))
        rows.append(_contig(f"e{i}", "TRB", str(i), f"CASE{i}", 6))
        rows.append(_contig(f"e{i}", "TRB", f"x{i}", f"CASSX{i}", 1))
    for i in range(20):
        rows.append(_contig(f"h{i}", "TRA", "9", "CAVPUB", 6))
        rows.append(_contig(f"h{i}", "TRB", f"h{i}", f"CASSH{i}", 6))
    for i in range(6):
        rows.append(_contig(f"p{i}", "TRA", "9", "CAVPUB", 1))
        rows.append(_contig(f"p{i}", "TRB", f"p{i}", f"CASSP{i}", 8))
    for i in range(8):
        rows.append(_contig(f"l{i}", "TRA", f"l{i}", f"CAL{i}", 8))
        rows.append(_contig(f"l{i}", "TRA", "9", "CAVPUB", 1))
        rows.append(_contig(f"l{i}", "TRB", f"l{i}", f"CASSL{i}", 8))
    endogenous = dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1))
    public = endogenous.filter(pl.col("cdr3") == "CAVPUB")
    primary_low = public.filter(pl.col("barcode").str.starts_with("p"))
    later = public.filter(pl.col("barcode").str.starts_with("l"))
    primary_high = public.filter(pl.col("barcode").str.starts_with("h"))
    assert primary_low.height == 6
    assert int(primary_low["endogenous"].sum()) == 6
    assert later.height == 8
    assert int(later["endogenous"].sum()) == 0
    assert int(primary_high["endogenous"].sum()) == 20


def test_skewed_receptor_cloud_keeps_a_mid_count():
    """A heavy tail is one negative binomial. A mid count stays endogenous.

    Two Poisson rates fit the same tail, and they would call the mid count a
    contaminant. The second rate has to beat the negative binomial.
    """
    rows = []
    for i in range(30):
        rows.append(_contig(f"e{i}", "TRA", str(i), f"CAE{i}", 6))
        rows.append(_contig(f"e{i}", "TRB", str(i), f"CASE{i}", 6))
        rows.append(_contig(f"e{i}", "TRB", f"x{i}", f"CASSX{i}", 1))
    primary = [1, 1, 1, 1, 4, 5, 5, 5, 5, 6, 7, 18, 26]
    for i, umi in enumerate(primary):
        rows.append(_contig(f"h{i}", "TRA", "9", "CAVSKEW", umi))
        rows.append(_contig(f"h{i}", "TRB", f"h{i}", f"CASSH{i}", 8))
    for i in range(8):
        rows.append(_contig(f"l{i}", "TRA", f"l{i}", f"CAL{i}", 8))
        rows.append(_contig(f"l{i}", "TRB", f"l{i}", f"CASSL{i}", 8))
        rows.append(_contig(f"l{i}", "TRB", "9", "CASSLEAK", 1))
    endogenous = dcode.mark_endogenous(dcode.chain_table(pl.DataFrame(rows), donor=1))
    skew = endogenous.filter(pl.col("cdr3") == "CAVSKEW")
    mid = skew.filter(pl.col("umi") == 5)
    high = skew.filter(pl.col("umi") == 26)
    assert mid.height == 4
    assert int(mid["endogenous"].sum()) == 4
    assert int(high["endogenous"].sum()) == 1


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


def _assert_posterior(score: float) -> None:
    assert 0.0 < score <= 1.0


def test_count_of_one_matches_background_and_eleven_is_specific():
    """Eleven UMIs is specific. One UMI is not, once one UMI is already common.

    A negative-control cell at 50 UMIs does not remove the call. The score is
    a posterior, and the hit is one cell.
    """
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    extras = [("donor_max", {ela: 0.0}, 50.0), ("ela", {ela: 11.0}, 0.0)]
    extras += [(f"one{i}", {ela: 1.0}, 0.0) for i in range(10)]
    clones, counts, controls = _singleton_panel(30, extras, [ela], nc_col)
    calls = dcode.call_specificities(clones, counts, [ela], controls, [nc_col])
    assert calls["clone_id"].to_list() == ["ela"]
    assert calls["peptide"].to_list() == ["ELAGIGILTV"]
    assert calls["mhc"].to_list() == ["HLA-A*02:01"]
    assert calls["n_cells"].to_list() == [1]
    _assert_posterior(float(calls["score"][0]))


def test_lone_count_above_a_zero_background_is_specific():
    """No UMI floor: the only cell above a zero background is the specific label."""
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        40, [("only", {ela: 1.0}, 0.0)], [ela], nc_col,
    )
    calls = dcode.call_specificities(clones, counts, [ela], controls, [nc_col])
    assert calls["clone_id"].to_list() == ["only"]
    assert calls["n_cells"].to_list() == [1]


def test_peptide_proportional_to_negative_control_is_background():
    """Repeating a background cell's peptide and control counts is not a call.

    Singletons sit on peptide = 2 x control. Four more cells at peptide 16 and
    control 8 match those singletons. Four cells at peptide 40 and control 0
    do not. The score is the specific posterior.
    """
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for level in (0.0, 1.0, 2.0, 4.0, 8.0):
        for i in range(8):
            barcode = f"b{int(level)}-{i}"
            count_rows.append({"barcode": barcode, pep: 2.0 * level})
            control_rows.append({"barcode": barcode, nc_col: level})
            clone_rows.append({"donor": 1, "barcode": barcode, "clone_id": barcode})
    for i in range(4):
        count_rows.append({"barcode": f"clean{i}", pep: 40.0})
        control_rows.append({"barcode": f"clean{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"clean{i}", "clone_id": "clean"})
        count_rows.append({"barcode": f"same{i}", pep: 16.0})
        control_rows.append({"barcode": f"same{i}", nc_col: 8.0})
        clone_rows.append({"donor": 1, "barcode": f"same{i}", "clone_id": "same"})
    calls = dcode.call_specificities(
        pl.DataFrame(clone_rows),
        pl.DataFrame(count_rows),
        [pep],
        pl.DataFrame(control_rows),
        [nc_col],
    )
    assert calls["clone_id"].to_list() == ["clean"]
    assert calls["n_cells"].to_list() == [4]
    _assert_posterior(float(calls["score"][0]))


def test_sticky_peptide_high_on_every_label_is_not_a_call():
    """KLGGALQAK at 20 on every label is the background mean, so it is not a call.

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
    _assert_posterior(float(calls["score"][0]))


def test_cells_at_the_background_rate_are_not_a_specific_clone():
    """Eight cells at 1 UMI match the singletons at 1. A cell at 30 does not.

    The eight cells are not summed into one count of 8. Sample size is the
    likelihood, and the score is not that sum.
    """
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(30):
        count_rows.append({"barcode": f"z{i}", pep: 1.0})
        control_rows.append({"barcode": f"z{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"z{i}", "clone_id": f"z{i}"})
    for i in range(10):
        count_rows.append({"barcode": f"o{i}", pep: 0.0})
        control_rows.append({"barcode": f"o{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"o{i}", "clone_id": f"o{i}"})
    for i in range(8):
        count_rows.append({"barcode": f"e{i}", pep: 1.0})
        control_rows.append({"barcode": f"e{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"e{i}", "clone_id": "eight"})
    count_rows.append({"barcode": "one", pep: 30.0})
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
    assert float(calls["score"][0]) != 8.0
    _assert_posterior(float(calls["score"][0]))


def test_median_does_not_assign_the_label():
    """(0, 100) has the higher log1p median and is one hot cell, not a call.

    (8, 8) is the specific label. A second label with the same counts gets the
    same posterior. Median and MAD are recorded on the check table afterwards.
    """
    pep = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nc_col = "A0201_ALIAPVHAV_NC"
    count_rows = []
    control_rows = []
    clone_rows = []
    for i in range(25):
        count_rows.append({"barcode": f"z{i}", pep: 0.0})
        control_rows.append({"barcode": f"z{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"z{i}", "clone_id": f"z{i}"})
    for i, value in enumerate((0.0, 100.0)):
        count_rows.append({"barcode": f"h{i}", pep: value})
        control_rows.append({"barcode": f"h{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"h{i}", "clone_id": "split"})
    for name in ("pair", "pair_copy"):
        for i in range(2):
            count_rows.append({"barcode": f"{name}{i}", pep: 8.0})
            control_rows.append({"barcode": f"{name}{i}", nc_col: 0.0})
            clone_rows.append({"donor": 1, "barcode": f"{name}{i}", "clone_id": name})
    clones = pl.DataFrame(clone_rows)
    counts = pl.DataFrame(count_rows)
    controls = pl.DataFrame(control_rows)
    split_median = float(np.median(np.log1p(np.array([0.0, 100.0]))))
    pair_median = float(np.log1p(8.0))
    assert split_median > pair_median
    calls, checks = dcode.assign_specificities(clones, counts, [pep], controls, [nc_col])
    assert set(calls["clone_id"].to_list()) == {"pair", "pair_copy"}
    scores = dict(zip(calls["clone_id"].to_list(), calls["score"].to_list(), strict=True))
    assert math.isclose(scores["pair"], scores["pair_copy"], rel_tol=0, abs_tol=1e-9)
    assert math.isclose(float(checks["median_log1p_called"][0]), pair_median, rel_tol=0, abs_tol=1e-9)
    assert float(checks["median_log1p_called"][0]) < split_median


def test_one_hot_cell_is_an_outlier_and_a_stable_clone_is_specific():
    """Four cells at 30 share a mean. One cell at 100 among three zeros does not."""
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
    _assert_posterior(float(calls["score"][0]))


def test_one_hot_label_is_not_specific_when_background_is_overdispersed():
    """Scattered counts of 1 and 2 keep the background overdispersed.

    Five zeros and one count of 80 are one hot cell. Four cells at 30 agree,
    and that label is the call. A varying negative control takes the scaled
    likelihood, and the same split holds.
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
    for i in range(15):
        count_rows.append({"barcode": f"a{i}", pep: 1.0})
        control_rows.append({"barcode": f"a{i}", nc_col: 0.0})
        clone_rows.append({"donor": 1, "barcode": f"a{i}", "clone_id": f"a{i}"})
    for i in range(8):
        count_rows.append({"barcode": f"c{i}", pep: 2.0})
        control_rows.append({"barcode": f"c{i}", nc_col: 1.0})
        clone_rows.append({"donor": 1, "barcode": f"c{i}", "clone_id": f"c{i}"})
    for i, value in enumerate((0.0, 0.0, 0.0, 0.0, 0.0, 80.0)):
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
    _assert_posterior(float(calls["score"][0]))


def test_cross_reactive_label_keeps_both_peptides():
    """A label elevated for two peptides keeps both. A count of 1 on a third does not."""
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    gil = "A0201_GILGFVFTL_Flu-MP_Influenza"
    nlv = "A0201_NLVPMVATV_pp65_CMV"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        30,
        [
            ("donor_max", {ela: 0.0, gil: 0.0, nlv: 0.0}, 50.0),
            ("both", {ela: 40.0, gil: 40.0, nlv: 1.0}, 0.0),
            *[(f"nlv{i}", {nlv: 1.0}, 0.0) for i in range(10)],
        ],
        [ela, gil, nlv],
        nc_col,
    )
    calls = dcode.call_specificities(clones, counts, [ela, gil, nlv], controls, [nc_col])
    assert set(calls["peptide"].to_list()) == {"ELAGIGILTV", "GILGFVFTL"}
    assert set(calls["clone_id"].to_list()) == {"both"}


def test_contig_umi_column_is_not_the_dextramer_score():
    """A contig-like UMI column of 500 is not read. The score is a posterior."""
    ela = "A0201_ELAGIGILTV_MART-1_Cancer"
    nc_col = "A0201_ALIAPVHAV_NC"
    clones, counts, controls = _singleton_panel(
        30,
        [("ela", {ela: 11.0}, 0.0), *[(f"one{i}", {ela: 1.0}, 0.0) for i in range(10)]],
        [ela],
        nc_col,
    )
    counts = counts.with_columns(pl.lit(500.0).alias("umis"))
    calls = dcode.call_specificities(clones, counts, [ela], controls, [nc_col])
    assert calls["clone_id"].to_list() == ["ela"]
    assert float(calls["score"][0]) != 500.0
    _assert_posterior(float(calls["score"][0]))


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
    """A leaked chain's contig UMI must not remove the beta-only dextramer call.

    The beta-only cells carry the binder beta plus a low contig-UMI copy of the
    paired clone's alpha. The score is the specific posterior, not log1p.
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
    assert float(hit["score"][0]) != math.log1p(100)
    _assert_posterior(float(hit["score"][0]))
