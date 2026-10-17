"""Curated parent labels and completed class-II pairs remain consistent across the built corpus."""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.emit.vdjdb3 import read_table

pytestmark = pytest.mark.release


def test_curated_metadata_is_complete_and_consistent():
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    if not (directory / "records.parquet").exists():
        pytest.skip("needs built records")
    records = read_table(directory, "records")
    missing = records.filter(
        pl.col("antigen.gene").is_null() | (pl.col("antigen.gene") == "")
        | pl.col("antigen.species").is_null() | (pl.col("antigen.species") == ""))
    assert missing.is_empty(), f"{missing.height} records need representative peptide parents"
    parents = records.select("antigen.epitope", "antigen.gene", "antigen.species").unique()
    conflicts = parents.group_by("antigen.epitope").len().filter(pl.col("len") > 1)
    assert conflicts.is_empty(), f"{conflicts.height} epitopes need one canonical parent pair"
    class_ii = records.filter(pl.col("mhc.class") == "MHCII")
    partial = class_ii.filter(
        pl.col("mhc.a").is_null() | (pl.col("mhc.a") == "")
        | pl.col("mhc.b").is_null() | (pl.col("mhc.b") == ""))
    assert partial.is_empty(), f"{partial.height} class-II records need compatible partner completion"
    hla = class_ii.filter(pl.col("mhc.a").str.starts_with("HLA-")
                         & pl.col("mhc.b").str.starts_with("HLA-"))
    wrong_family = hla.filter(
        pl.col("mhc.a").str.extract(r"HLA-(DP|DQ|DR)", 1)
        != pl.col("mhc.b").str.extract(r"HLA-(DP|DQ|DR)", 1))
    assert wrong_family.is_empty(), f"{wrong_family.height} class-II records mix incompatible loci"
    wrong_role = hla.filter(
        ~pl.col("mhc.a").str.contains(r"^HLA-(DRA|DQA[12]|DPA[12])\*")
        | ~pl.col("mhc.b").str.contains(r"^HLA-(DRB[1345]|DQB[12]|DPB[12])\*"))
    assert wrong_role.is_empty(), f"{wrong_role.height} class-II records have incorrect chain roles"
