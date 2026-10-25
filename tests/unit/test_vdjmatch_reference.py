"""The published consumer must read the user-selected build without downloading a release."""

import polars as pl
import pytest
from vdjmatch.db import vdjdb


def test_local_build_reference_is_current_and_offline(tmp_path, monkeypatch):
    def no_download(**kwargs):
        pytest.fail("a local build must not fetch a remote reference")

    monkeypatch.setattr(vdjdb, "fetch_latest", no_download)
    tables = tmp_path / "chosen-build" / "tables"
    tables.mkdir(parents=True)
    record = {
        "record_id": "VDJDB0000000001", "gene": "TRB", "cdr3": "CASSIRSSYEQYF",
        "v_segm": "TRBV19*01", "j_segm": "TRBJ2-7*01", "species": "HomoSapiens",
        "antigen_epitope": "GILGFVFTL", "mhc_a": "HLA-A*02:01", "mhc_b": "B2M",
        "mhc_class": "MHCI", "reference_id": "PMID:1", "vdjdb_score": 1,
    }
    path = tables / "vdjdb.parquet"
    pl.DataFrame([record]).write_parquet(path)
    assert vdjdb.load(tables)["record_id"].to_list() == [record["record_id"]]
    record["record_id"] = "VDJDB0000000002"
    pl.DataFrame([record]).write_parquet(path)
    assert vdjdb.load(tables)["record_id"].to_list() == [record["record_id"]]
