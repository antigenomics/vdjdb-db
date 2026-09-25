"""The proprietary-data guard. TCRvdb must never reach the repository or a release."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from vdjdb.validate import scan_paths, tcrvdb_path
from vdjdb.validate.guard import ENV_VAR


def test_clean_tree_has_no_leaks(tmp_path):
    (tmp_path / "vdjdb.txt").write_text("gene\tcdr3\nTRB\tCASSF\n")
    assert scan_paths(list(tmp_path.iterdir()), root=tmp_path) == []


@pytest.mark.parametrize("name", [
    "01_05_2025_TCRvdb.csv",
    "tcrvdb_labels.tsv",
    "matchmakers.parquet",
    "MatchMakers_export.csv",
])
def test_filename_fingerprints_are_caught(tmp_path, name):
    (tmp_path / name).write_text("a,b\n1,2\n")
    leaks = scan_paths(list(tmp_path.iterdir()), root=tmp_path)
    assert len(leaks) == 1 and "filename" in leaks[0].reason


def test_a_renamed_copy_is_caught_by_its_header(tmp_path):
    """Renaming the file must not defeat the guard."""
    p = tmp_path / "harmless_numbers.csv"
    p.write_text("name,clonotype_aa,epitope_aa,hla_short,baseMean,padj\nx,CASSF,GILGFVFTL,A2,1.0,0.01\n")
    leaks = scan_paths([p], root=tmp_path)
    assert len(leaks) == 1 and "fingerprint" in leaks[0].reason


def test_a_file_with_only_some_fingerprint_columns_is_not_flagged(tmp_path):
    p = tmp_path / "ours.tsv"
    p.write_text("epitope_aa\tpadj\nGILGFVFTL\t0.01\n")
    assert scan_paths([p], root=tmp_path) == []


def test_tcrvdb_path_refuses_to_guess(monkeypatch):
    monkeypatch.delenv(ENV_VAR, raising=False)
    with pytest.raises(RuntimeError, match="proprietary"):
        tcrvdb_path()


def test_tcrvdb_path_reports_a_missing_file(monkeypatch, tmp_path):
    monkeypatch.setenv(ENV_VAR, str(tmp_path / "nope.csv"))
    with pytest.raises(FileNotFoundError):
        tcrvdb_path()


def test_tcrvdb_path_returns_the_configured_file(monkeypatch, tmp_path):
    p = tmp_path / "held_out.csv"
    p.write_text("x\n")
    monkeypatch.setenv(ENV_VAR, str(p))
    assert tcrvdb_path() == p


def test_the_real_repository_is_clean():
    """The guard CI runs. If this fails, something proprietary has been committed."""
    root = Path(__file__).resolve().parents[2]
    tracked = [
        root / line
        for line in os.popen(f"git -C {root} ls-files").read().splitlines()
        if line
    ]
    leaks = scan_paths([p for p in tracked if p.exists()], root=root)
    assert leaks == [], "proprietary data in the repository: " + "; ".join(str(x) for x in leaks)
