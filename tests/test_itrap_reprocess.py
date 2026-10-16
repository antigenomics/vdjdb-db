"""Reshape deposited ITRAP counts into the committed caller. No cutoff grid."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import polars as pl

_SCRIPT = Path(__file__).resolve().parents[1] / "processing" / "itrap_reprocess.py"
_spec = importlib.util.spec_from_file_location("itrap_reprocess", _SCRIPT)
assert _spec is not None and _spec.loader is not None
itrap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(itrap)


def _gem(gem: str, donor: str, alpha: str, beta: str, peptides: str, umis: str) -> dict:
    return {
        "gem": gem,
        "sample_id": donor,
        "genes_lst_TRA": "TRAV12-2;TRAJ1;TRAC",
        "cdr3_lst_TRA": alpha,
        "umi_count_lst_TRA": "5",
        "genes_lst_TRB": "TRBV19;TRBD1;TRBJ1-1;TRBC1",
        "cdr3_lst_TRB": beta,
        "umi_count_lst_TRB": "6",
        "peptide_HLA_lst": peptides,
        "umi_count_lst_mhc": umis,
    }


def test_mismatched_chain_list_is_skipped_and_peptides_stay():
    frame = pl.DataFrame(
        [
            _gem(
                "g1",
                "2.0",
                "CAVNNARLMF",
                "CASSIRSSYEQYF",
                "GILGFVFTL A0201|NLVPMVATV A0201",
                "4.0|1.0",
            ),
        ]
    )
    frame = frame.with_columns(pl.lit("CASSIRSSYEQYF|CASSSTOOOSHORT").alias("cdr3_lst_TRB"))
    frame = frame.with_columns(pl.lit("6").alias("umi_count_lst_TRB"))
    contigs, features, stats = itrap.itrap1_tables(frame)
    assert stats["skipped_chain_lists"] == 1
    assert contigs.filter(pl.col("chain") == "TRB").height == 0
    assert "A0201_GILGFVFTL" in features.columns
    assert features["A0201_GILGFVFTL"][0] == 4.0


def test_two_clones_keep_their_own_peptide():
    rows = []
    for index, gil, nlv in ((1, "30", "0"), (2, "28", "0"), (3, "32", "1")):
        rows.append(
            _gem(
                f"a{index}",
                "1.0",
                "CAVNNARLMF",
                "CASSIRSSYEQYF",
                "GILGFVFTL A0201|NLVPMVATV A0201",
                f"{gil}|{nlv}",
            )
        )
    for index, gil, nlv in ((1, "0", "25"), (2, "1", "27"), (3, "0", "22")):
        rows.append(
            _gem(
                f"b{index}",
                "1.0",
                "CAVGVRGGGNKLTF",
                "CASSVTGPPREDGYTF",
                "GILGFVFTL A0201|NLVPMVATV A0201",
                f"{gil}|{nlv}",
            )
        )
    contigs, features, _stats = itrap.itrap1_tables(pl.DataFrame(rows))
    result = itrap.run_frames(contigs, features)
    called = {(row[0], row[1]) for row in result["calls"].select("clone_id", "peptide").iter_rows()}
    gil_clone = [key for key in called if key[1] == "GILGFVFTL"]
    nlv_clone = [key for key in called if key[1] == "NLVPMVATV"]
    assert len(gil_clone) == 1
    assert len(nlv_clone) == 1
    assert gil_clone[0][0] != nlv_clone[0][0]
