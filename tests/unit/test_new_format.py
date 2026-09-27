"""The new format: the joined view, the parquet round-trip, and the generated schema."""
from __future__ import annotations

import json

import polars as pl
import pytest

from vdjdb.assemble.evidence import build_evidence
from vdjdb.emit import legacy, vdjdb3
from vdjdb.schema import CHAIN_COLUMNS, RECORD_COLUMNS, TABLES

#: Columns whose dtype the legacy projection depends on: the rest are strings.
_CHAIN_DTYPES = {"v.end": pl.Int64, "j.start": pl.Int64,
                 "fix.needed": pl.Boolean, "fix.good": pl.Boolean,
                 "v.canonical": pl.Boolean, "j.canonical": pl.Boolean,
                 "cdr3nt.pgen": pl.Float64, "cdr3nt.margin": pl.Float64}

RECORDS = [
    {"record_id": "VDJDB0000000001", "species": "HomoSapiens", "mhc.a": "HLA-A*02:01",
     "mhc.b": "B2M", "mhc.class": "MHCI", "antigen.epitope": "GILGFVFTL", "antigen.gene": "M",
     "antigen.species": "InfluenzaA", "reference.id": "PMID:1", "vdjdb.score": 2,
     "chunk.file": "PMID_1.txt", "chunk.row": 0, "method.identification": "tetramer-sort"},
    {"record_id": "VDJDB0000000002", "species": "HomoSapiens", "mhc.a": "HLA-A*02:01",
     "mhc.b": "B2M", "mhc.class": "MHCI", "antigen.epitope": "GILGFVFTL", "antigen.gene": "M",
     "antigen.species": "InfluenzaA", "reference.id": "PMID:2", "vdjdb.score": 1,
     "chunk.file": "PMID_2.txt", "chunk.row": 0, "method.identification": "culture"},
]

CHAINS = [
    {"record_id": "VDJDB0000000001", "gene": "TRB", "clonotype_id": "CT2df51ea0980263d2",
     "cdr3": "CASSIRSSYEQYF",
     "v.segm": "TRBV10-3*01", "j.segm": "TRBJ2-7*01", "v.end": 4, "j.start": 8,
     "cdr3.original": "CASSIRSSYEQYF", "TCR_hash": "abc"},
    {"record_id": "VDJDB0000000002", "gene": "TRB", "clonotype_id": "CT2df51ea0980263d2",
     "cdr3": "CASSIRSSYEQYF",
     "v.segm": "TRBV10-3*01", "j.segm": "TRBJ2-7*01", "v.end": 4, "j.start": 8,
     "cdr3.original": "CASSIRSSYEQYF", "TCR_hash": "abc"},
]


def _frame(rows: list[dict], columns: tuple[str, ...], dtypes: dict) -> pl.DataFrame:
    blank = {pl.Int64: 0, pl.Boolean: False, pl.Float64: None}
    filled = [{c: r.get(c, blank.get(dtypes.get(c), "")) for c in columns} for r in rows]
    return pl.DataFrame(filled, schema={c: dtypes.get(c, pl.String) for c in columns})


@pytest.fixture
def tables() -> dict[str, pl.DataFrame]:
    from vdjdb.assemble.epitopes import build_epitopes, build_restriction

    records = _frame(RECORDS, RECORD_COLUMNS, {"vdjdb.score": pl.Int64, "chunk.row": pl.Int64})
    chains = _frame(CHAINS, CHAIN_COLUMNS, _CHAIN_DTYPES)
    return {"records": records, "chains": chains,
            "evidence": build_evidence(records, chains, release="v1"),
            "epitopes": build_epitopes(records, chains),
            "restriction": build_restriction(records)}


def test_the_view_declares_every_evidence_column_even_without_a_producer(tables):
    """The view's shape must not change as producers land -- absent evidence is `false`, not absent."""
    view = vdjdb3.joined(tables)
    for col in vdjdb3.VIEW_EVIDENCE_COLUMNS:
        assert view.schema[col] == pl.Boolean, col
    assert view["evidence.motif.tcrnet"].to_list() == [False, False]
    assert view["evidence.validation.same.study"].to_list() == [False, False]


def test_the_view_is_one_row_per_chain_and_carries_the_evidence_it_has(tables):
    view = vdjdb3.joined(tables)
    assert view.height == tables["chains"].height
    # two papers report this clonotype against this epitope, so both records carry the evidence
    assert view["evidence.validation.independent"].to_list() == [True, True]


def test_record_level_evidence_is_refused_rather_than_silently_dropped(tables):
    """Structure evidence covers a receptor, not a chain. When it lands it needs a second join."""
    tables["evidence"] = tables["evidence"].with_columns(pl.lit("").alias("gene"))
    with pytest.raises(NotImplementedError, match="record_id"):
        vdjdb3.joined(tables)


def test_tables_round_trip_through_parquet_unchanged(tmp_path, tables):
    vdjdb3.write_all(tables, tmp_path)
    back = vdjdb3.read_tables(tmp_path)
    for name in ("records", "chains", "evidence"):
        assert back[name].equals(tables[name]), name


def test_the_legacy_export_is_a_projection_of_what_shipped(tmp_path, tables):
    """Phase 6's closing criterion in miniature: `make legacy` reads the tables, never `chunks/`.

    On the real corpus the three projected files are byte-identical to the in-memory build and the
    release comparison passes; here the same property is asserted cheaply on every run.
    """
    direct = legacy.write_all(tables, tmp_path / "direct")
    vdjdb3.write_all(tables, tmp_path / "new")
    projected = legacy.write_all(vdjdb3.read_tables(tmp_path / "new"), tmp_path / "projected")
    assert set(direct) == set(projected)
    for name in direct:
        assert direct[name].read_bytes() == projected[name].read_bytes(), name


def test_the_generated_schema_describes_every_column_of_every_table(tmp_path, tables):
    vdjdb3.write_all(tables, tmp_path)
    doc = json.loads((tmp_path / "vdjdb.schema.json").read_text())
    described = {f["name"]: f for f in doc["fields"]}
    for table, columns in TABLES.items():
        assert doc["tables"][table] == list(columns)
        for i, col in enumerate(columns):
            assert described[col]["position"][table] == i, (table, col)
    # dtypes are read off the written frames, so the schema cannot claim a type the files lack
    assert described["record_id"]["dtype"] == "String"
    assert described["evidence_score"]["dtype"] == "Float64"
    # A string, not a uint: the id is `CT` plus 16 hex digits of a sha256 we own, because polars
    # does not specify `Expr.hash` across versions and this id ships (ROADMAP.md section 10.3).
    assert described["clonotype_id"]["dtype"] == "String"
