"""Reshape deposited ITRAP counts into the committed caller. No cutoff grid."""

from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import polars as pl
import pytest

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
    assert set(result) == {"chains", "edges", "clones", "calls", "checks"}


def test_write_result_uses_the_caller_table_names(tmp_path: Path) -> None:
    """The ITRAP reshaper writes the same five tables the 10x command writes."""
    result = {
        key: pl.DataFrame({"donor": [1], "barcode": ["g1"]})
        for key in ("chains", "edges", "clones", "calls", "checks")
    }
    itrap.write_result(result, tmp_path)
    written = {path.name for path in tmp_path.iterdir()}
    assert written == {"chains.tsv", "edges.tsv", "clones.tsv", "calls.tsv", "checks.tsv"}
    for path in tmp_path.iterdir():
        assert pl.read_csv(path, separator="\t").columns == ["donor", "barcode"]


def test_reused_barcodes_keep_feature_counts_with_their_donor():
    rows = []
    for donor in ("1.0", "2.0"):
        for index in range(3):
            rows.append(_gem(
                f"g{index}", donor, "CAVNNARLMF", "CASSIRSSYEQYF",
                "GILGFVFTL A0201", "30" if donor == "1.0" else "0",
            ))
    contigs, features, _ = itrap.itrap1_tables(pl.DataFrame(rows))
    result = itrap.run_frames(contigs, features)
    checks = result["checks"].filter(pl.col("donor") == 2)
    assert checks["mu0"].to_list() == [0.0]
    assert checks["n_called"].to_list() == [0]


def test_repeated_runs_write_identical_bytes(tmp_path: Path):
    rows = [_gem(f"g{i}", "1.0", f"CAV{aa}NARLMF", "CASSIRSSYEQYF", "GILGFVFTL A0201", "0")
            for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")]
    contigs, features, _ = itrap.itrap1_tables(pl.DataFrame(rows))
    previous = None
    for _ in range(3):
        itrap.write_result(itrap.run_frames(contigs, features), tmp_path)
        current = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
        if previous is not None:
            assert current == previous
        previous = current


def _matrix_files(tmp_path: Path, body: str):
    paths = [tmp_path / name for name in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")]
    for path, text in zip(paths, ("%%MatrixMarket matrix coordinate integer general\n" + body,
                                 "BC-1\n", "probe\tprobe\tAntibody Capture\n"), strict=True):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    return paths


@pytest.mark.parametrize("body", [
    "1 1 1\n1 0 3\n", "1 1 1\n0 1 3\n", "1 1 1\n2 1 3\n",
    "1 2 1\n1 1 3\n", "1 1 2\n1 1 3\n", "1 1 1\n1 1 0.5\n",
])
def test_invalid_matrix_market_input_is_rejected(tmp_path: Path, body: str):
    paths = _matrix_files(tmp_path, body)
    with pytest.raises(ValueError, match="MatrixMarket"):
        itrap._matrix_antibody(*paths, {"probe": "A0201_GILGFVFTL"})


def test_matrix_market_entries_sum_for_the_mapped_peptide(tmp_path: Path):
    paths = _matrix_files(tmp_path, "1 1 2\n1 1 3\n1 1 4\n")
    frame = itrap._matrix_antibody(*paths, {"probe": "A0201_GILGFVFTL"})
    assert frame.to_dicts() == [{"barcode": "BC-1", "A0201_GILGFVFTL": 7.0}]


def test_viral_loader_with_no_contigs_returns_typed_empty_frames(tmp_path: Path):
    matrix, barcodes, features = _matrix_files(tmp_path, "1 1 0\n")
    alpha, beta, meta, key = [tmp_path / name for name in ("alpha.tsv", "beta.tsv", "meta.tsv", "key.tsv")]
    for path in (alpha, beta):
        path.write_text("cell_id\tproductive\tv_call\tj_call\tjunction_aa\tumi_count\n")
    meta.write_text("cell\tBC\n1OS_BC\tdonor1\n")
    key.write_text("Barcode\tSequence\tHLA\nprobe\tGILGFVFTL\tA0201\n")
    contigs, counts, stats = itrap.viral_tables(alpha, beta, meta, matrix, barcodes, features, key)
    assert contigs.is_empty() and counts.is_empty()
    assert contigs.schema["donor"] == pl.Int64
    assert stats["contigs"] == stats["donors"] == 0
