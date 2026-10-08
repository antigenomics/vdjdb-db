"""The assessment covers reported observations without changing their identity."""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.assemble.assessment import KEY, PROVENANCE, specification
from vdjdb.emit.vdjdb3 import read_table
from vdjdb.schema import EPITOPE_ASSESSMENT_COLUMNS

pytestmark = pytest.mark.release


@pytest.fixture(scope="module")
def frames():
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    if not (directory / "epitope_assessment.parquet").exists():
        assert not os.environ.get("CI"), "CI must build the epitope assessment table"
        pytest.skip("no built epitope assessment table")
    return read_table(directory, "records"), read_table(directory, "epitope_assessment")


def test_reported_pairs_preserve_provenance_and_support(frames):
    records, assessment = frames
    source_key = [c for c in KEY if c != "mhc.species"]
    expected = (records.group_by(source_key)
                .agg(pl.len().alias("records"), pl.col("reference.id").n_unique().alias("references"))
                .sort(source_key))
    reported = assessment.filter(pl.col("reported")).select(*source_key, "records", "references")
    assert reported.sort(source_key).equals(expected)
    assert assessment.select(*KEY, "prediction.allele").unique().height == assessment.height
    assert tuple(assessment.columns) == EPITOPE_ASSESSMENT_COLUMNS
    assert sum(assessment.null_count().row(0)) == 0
    inferred = assessment.filter(~pl.col("reported"))
    assert inferred.filter((pl.col("mhc.a") != "") | (pl.col("mhc.b") != "") |
                           (pl.col("records") != 0) | (pl.col("references") != 0)).is_empty()
    assert inferred.filter(~pl.col("prediction.best") &
                           (pl.col("presentation.band") == "non-binder")).is_empty()


def test_scored_core_coordinates_and_model_provenance(frames):
    _, assessment = frames
    if "reference_not_supplied" in assessment["assessment.status"]:
        assert not os.environ.get("CI"), "CI must supply the pinned presentation reference"
        pytest.skip("no presentation reference supplied to this local build")
    scored = assessment.filter(pl.col("assessment.status") == "scored")
    assert not scored.is_empty()
    spec = specification()
    assert scored["reference.sha256"].unique().to_list() == [spec["sha256"]]
    assert scored["mhcmatch.version"].unique().to_list() == [spec["mhcmatch_version"]]
    assert scored["prediction.background"].unique().to_list() == [spec["background"]]
    assert scored["prediction.footprint"].unique().to_list() == [spec["footprint"]]
    assert scored.filter(pl.col("core").str.len_chars()
                         != pl.col("core.tcr.facing").str.len_chars()).is_empty()
    two = scored.filter(pl.col("mhc.class") == "MHCII")
    assert two.filter(pl.col("core") != pl.col("prediction.peptide").str.slice(
        pl.col("core.offset").cast(pl.Int64), 9)).is_empty()
    offset = pl.col("prediction.offset").cast(pl.Int64)
    assert scored.filter(pl.col("prediction.peptide") != pl.col("antigen.epitope").str.slice(
        offset, pl.col("prediction.peptide").str.len_chars())).is_empty()
    best = scored.filter(pl.col("prediction.best")).select(*PROVENANCE, "prediction.allele").unique()
    assert best.group_by(PROVENANCE).len().filter(pl.col("len") != 1).is_empty()
