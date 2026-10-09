"""Reported observations survive predictions, flank comparisons and missing coverage."""
from __future__ import annotations

import hashlib
import os

import polars as pl
import pytest

from vdjdb.assemble import assessment as A
from vdjdb.schema import EPITOPE_ASSESSMENT_COLUMNS


def records() -> pl.DataFrame:
    return pl.DataFrame({
        "antigen.epitope": ["PKYVKQNTLKLAT", "PKYVKQNTLKLAT", "PKYVKQNTLKLAT", "AA"],
        "antigen.species": ["InfluenzaA"] * 4,
        "antigen.gene": ["HA", "HA", "Other", "HA"],
        "species": ["HomoSapiens"] * 4, "mhc.class": ["MHCII"] * 4,
        "mhc.a": ["HLA-DRA*01:01"] * 4,
        "mhc.b": ["HLA-DRB1*01:01", "HLA-DRB1*01:01", "HLA-DRB1*04:01", "HLA-DRB1*01:01"],
        "reference.id": ["PMID:1", "PMID:2", "PMID:3", "PMID:4"],
    })


def test_no_reference_retains_each_pair_and_parent_gene():
    original = records()
    result = A.build_assessment(original)
    assert result.height == 3
    assert tuple(result.columns) == EPITOPE_ASSESSMENT_COLUMNS
    assert result["reported"].all()
    assert result["records"].sum() == original.height
    assert result["assessment.status"].unique().to_list() == ["reference_not_supplied"]
    assert sum(result.null_count().row(0)) == 0
    assert result.equals(A.build_assessment(original.reverse()))


def test_empty_input_and_invalid_jobs():
    assert A.build_assessment(records().head(0)).height == 0
    with pytest.raises(ValueError, match="positive"):
        A.build_assessment(records(), jobs=0)


def test_checksum_failure_precedes_any_scoring(tmp_path):
    bad = tmp_path / "reference.tsv"
    bad.write_text("not the declared reference")
    with pytest.raises(ValueError, match="checksum"):
        A.build_assessment(records(), reference=bad)


@pytest.mark.parametrize("band, expected", [("weak", 1), ("non-binder", 0)])
def test_join_keeps_predictions_for_other_provenance_and_never_invents_support(
        monkeypatch, tmp_path, band, expected):
    def score(task):
        return [
            {"antigen.epitope": "PKYVKQNTLKLAT", "mhc.a": "HLA-DRA*01:01",
             "mhc.b": "HLA-DRB1*01:01", "reported": True, "prediction.allele": "DRB1_0101",
             "core": "YVKQNTLKL", "core.offset": "2", "assessment.status": "scored",
             "presentation.band": band},
            {"antigen.epitope": "PKYVKQNTLKLAT", "reported": False,
             "prediction.allele": "DRB1_0101", "core": "YVKQNTLKL", "core.offset": "2",
             "assessment.status": "scored", "presentation.band": band},
            {"antigen.epitope": "AA", "mhc.a": "HLA-DRA*01:01", "mhc.b": "HLA-DRB1*01:01",
             "reported": True, "assessment.status": "unsupported_peptide"},
        ]
    monkeypatch.setattr(A, "verify_reference", lambda *_: None)
    monkeypatch.setattr(A, "_score_group", score)
    result = A.build_assessment(records(), reference=tmp_path / "reference")
    reported = result.filter(pl.col("reported"))
    assert reported["records"].sum() == 4
    assert reported.filter(pl.col("antigen.gene") == "Other")["mhc.b"].to_list() == ["HLA-DRB1*04:01"]
    predicted = result.filter(~pl.col("reported"))
    assert predicted.height == expected
    assert predicted["antigen.gene"].to_list() == ["Other"] * expected
    assert predicted["records"].to_list() == [0] * expected
    assert predicted["references"].to_list() == [0] * expected
    assert sum(result.null_count().row(0)) == 0


def test_environment_restored_after_failure(monkeypatch):
    monkeypatch.setenv("MHCMATCH_CALIBRATION_CACHE", "before")
    monkeypatch.setenv("OMP_NUM_THREADS", "7")
    with pytest.raises(RuntimeError), A._environment(parallel=True):
        assert os.environ["MHCMATCH_CALIBRATION_CACHE"] == "off"
        assert os.environ["OMP_NUM_THREADS"] == "1"
        raise RuntimeError("failed scorer")
    assert os.environ["MHCMATCH_CALIBRATION_CACHE"] == "before"
    assert os.environ["OMP_NUM_THREADS"] == "7"


def test_class_two_uses_allele_register_and_preserves_incomplete_dq(monkeypatch):
    import mhcmatch
    import mhcmatch.predict

    class Model:
        def score(self, peptide, allele):
            return 1.0

        def best_register(self, peptide, allele):
            return 3, 1.0

    class Calibration:
        def percent_rank(self, allele, score, length=None):
            assert length in (13, 14)
            return 1.0

        def p_present(self, allele, score):
            return 0.75

    real = mhcmatch.Store.from_records([
        {"epitope": "PKYVKQNTLKLAT", "mhc_class": "II", "mhc_a": "HLA-DRA*01:01",
         "mhc_b": "HLA-DRB1*01:01"},
    ])
    monkeypatch.setattr(mhcmatch.Store, "from_pmhc", lambda **_: real)
    monkeypatch.setattr(mhcmatch.predict, "build_scorer", lambda *_, **__: (Model(), Calibration(), None))
    got = A._score_group(("unused", A.specification(), "human", "mhc2", [
        ("PKYVKQNTLKLAT", "HLA-DRA*01:01", "HLA-DRB1*01:01"),
        ("PKYVKQNTLKLAT", "", "HLA-DQB1*03:01"),
        ("AA", "HLA-DRA*01:01", "HLA-DRB1*01:01"),
        ("PKYVKQNTLKLATA", "HLA-DRA*01:01", "HLA-DRB1*01:01"),
    ]))
    reported = next(row for row in got if row["mhc.b"] == "HLA-DRB1*01:01"
                    and row["assessment.status"] == "scored")
    assert reported["core"] == "VKQNTLKLA"
    assert reported["core.offset"] == "3"
    assert reported["core.tcr.facing"] == "XKQXTXKLX"  # mhcmatch class-default anchors
    partial = next(row for row in got if row["mhc.b"] == "HLA-DQB1*03:01")
    assert partial["assessment.status"] == "allele_not_in_panel"
    assert "core" not in partial
    unsupported = next(row for row in got if row["antigen.epitope"] == "AA")
    assert unsupported["assessment.status"] == "unsupported_peptide"
    assert unsupported["prediction.allele"] == reported["prediction.allele"]
    assert unsupported["allele.resolution"] == "exact"

    monkeypatch.setattr(Model, "score", lambda *args: float("nan"))
    unscored = A._score_group(("unused", A.specification(), "human", "mhc2", [
        ("PKYVKQNTLKLAT", "HLA-DRA*01:01", "HLA-DRB1*01:01"),
    ]))[0]
    assert unscored["assessment.status"] == "not_scorable"
    assert unscored["allele.resolution"] == "exact"
    assert unscored["prediction.allele"] == reported["prediction.allele"]


def test_class_one_footprint_is_not_a_contiguous_subsequence():
    from mhcmatch.store import binding_core

    assert binding_core("GILGFVFTLA", "mhc1") == ("GILGFFTLA", 0)


def test_mhc_species_is_independent_of_receptor_species():
    got = A.build_assessment(records().with_columns(pl.lit("MusMusculus").alias("species")))
    assert got["species"].unique().to_list() == ["MusMusculus"]
    assert got["mhc.species"].unique().to_list() == ["HomoSapiens"]


def test_real_scorer_serial_parallel_and_input_order_are_identical(tmp_path, monkeypatch):
    # A small reference fixture exercises the published scorer and process initialization offline.
    reference = tmp_path / "pmhc.tsv"
    reference.write_text(
        "epitope\tmhc_a\tmhc_b\tmhc_class\tmhc_species\n"
        "GILGFVFTL\tHLA-A*02:01\tB2M\tI\tHomoSapiens\n"
        "NLVPMVATV\tHLA-A*02:01\tB2M\tI\tHomoSapiens\n"
        "PKYVKQNTLKLAT\tHLA-DRA*01:01\tHLA-DRB1*01:01\tII\tHomoSapiens\n"
        "YVKQNTLKLAT\tHLA-DRA*01:01\tHLA-DRB1*01:01\tII\tHomoSapiens\n"
    )
    spec = {**A.specification(), "sha256": hashlib.sha256(reference.read_bytes()).hexdigest()}
    monkeypatch.setattr(A, "specification", lambda: spec)
    first = records().head(1)
    second = first.with_columns(pl.lit("GILGFVFTL").alias("antigen.epitope"),
                                pl.lit("MHCI").alias("mhc.class"),
                                pl.lit("HLA-A*02:01").alias("mhc.a"),
                                pl.lit("B2M").alias("mhc.b"))
    inputs = pl.concat([first, second])
    serial = A.build_assessment(inputs, reference=reference, jobs=1)
    parallel = A.build_assessment(inputs.reverse(), reference=reference, jobs=2)
    assert serial.equals(parallel)
    assert serial.filter(pl.col("reported"))["assessment.status"].to_list() == ["scored", "scored"]
    assert serial.write_csv(separator="\t") == parallel.write_csv(separator="\t")
