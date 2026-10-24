"""The new format: the joined view, the parquet round-trip, and the generated schema."""
from __future__ import annotations

import gzip
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
                 "cdr3nt.pgen": pl.Float64, "cdr3nt.margin": pl.Float64,
                 **{f"{b}.{u}": pl.Int64 for b in ("v.end", "j.start", "d.start", "d.end")
                    for u in ("nt", "aa")}}

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
    from vdjdb.assemble.assessment import build_assessment
    from vdjdb.assemble.epitopes import build_epitopes, build_restriction

    records = _frame(RECORDS, RECORD_COLUMNS, {"vdjdb.score": pl.Int64, "chunk.row": pl.Int64})
    chains = _frame(CHAINS, CHAIN_COLUMNS, _CHAIN_DTYPES)
    return {"records": records, "chains": chains,
            "evidence": build_evidence(records, chains, release="v1"),
            "epitopes": build_epitopes(records, chains),
            "restriction": build_restriction(records), "epitope_assessment": build_assessment(records)}


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


def test_record_level_evidence_reaches_both_chains_and_zero_motifs_stay_false(tables):
    tables['evidence'] = tables['evidence'].with_columns(
        pl.lit('').alias('gene'), pl.lit('structure_model').alias('evidence_type'))
    assert vdjdb3.joined(tables)['evidence.structure.model'].to_list() == [True, True]
    tables['evidence'] = tables['evidence'].with_columns(
        pl.lit('TRB').alias('gene'), pl.lit('motif_tcrnet').alias('evidence_type'),
        pl.lit(0.).alias('evidence_score'))
    assert vdjdb3.joined(tables)['evidence.motif.tcrnet'].to_list() == [False, False]


def test_tables_round_trip_through_parquet_unchanged(tmp_path, tables):
    vdjdb3.write_all(tables, tmp_path)
    back = vdjdb3.read_tables(tmp_path)
    for name in tables:
        assert back[name].equals(tables[name]), name


def test_older_tables_without_assessment_still_project_to_legacy(tmp_path, tables):
    vdjdb3.write_all(tables, tmp_path / "tables")
    for suffix in ("parquet", "tsv.gz"):
        (tmp_path / "tables" / f"epitope_assessment.{suffix}").unlink()
    old = vdjdb3.read_tables(tmp_path / "tables")
    assert "epitope_assessment" not in old
    legacy.write_all(old, tmp_path / "legacy")


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


def test_the_tidy_tables_ship_underscore_names_and_the_legacy_ones_do_not(tmp_path, tables):
    """The convention, asserted on the files rather than on the mapping.

    A dot is a table qualifier in SQL and blocks attribute access in most dataframe libraries. The
    legacy tables cannot follow: their column order is a positional contract `vdjdb-web` parses, and
    `vdjdb.meta.txt` must match `vdjdb.txt` column for column.
    """
    written = vdjdb3.write_all(tables, tmp_path / "new")
    for name in written:
        if not name.endswith(".tsv.gz"):
            continue
        header = gzip.decompress((tmp_path / "new" / name).read_bytes()).decode().splitlines()[0].split("\t")
        assert not [c for c in header if "." in c], (name, [c for c in header if "." in c])

    legacy_files = legacy.write_all(tables, tmp_path / "legacy")
    head = legacy_files["vdjdb.txt"].read_text().split("\n", 1)[0].split("\t")
    assert "antigen.epitope" in head and "antigen_epitope" not in head


def test_reading_a_tidy_table_back_restores_every_internal_name(tmp_path, tables):
    """The rename has to be exactly invertible or a round trip drops a column silently.

    `write_all` is the only place the tidy names are applied and `read_table` the only place they
    are undone, so this is the property that lets everything inside `assemble/` keep the dotted
    names.
    """
    vdjdb3.write_all(tables, tmp_path)
    for name, frame in tables.items():
        back = vdjdb3.read_table(tmp_path, name)
        assert back.columns == frame.select(back.columns).columns
        assert set(back.columns) == set(frame.columns), name

    # and a projection asks in internal names, not shipped ones
    got = vdjdb3.read_table(tmp_path, "records", ["record_id", "antigen.epitope"])
    assert got.columns == ["record_id", "antigen.epitope"]


def test_compressed_text_is_deterministic_and_matches_parquet(tmp_path, tables):
    import io
    for folder in ('a', 'b'):
        vdjdb3.write_all(tables, tmp_path/folder)
    for name in tables:
        a = tmp_path/'a'/f'{name}.tsv.gz'
        assert a.read_bytes() == (tmp_path/'b'/a.name).read_bytes()
        text = pl.read_csv(io.BytesIO(gzip.decompress(a.read_bytes())), separator='\t')
        assert text.height == tables[name].height
    assert not (tmp_path/'a'/'vdjdb.tsv.gz').exists()


def test_motif_run_replaces_its_evidence_and_refreshes_view(tmp_path, tables, monkeypatch):
    from vdjdb import motifs
    from vdjdb.assemble.evidence import empty
    vdjdb3.write_all(tables, tmp_path/'tables')
    ev = empty().vstack(pl.DataFrame({
        'record_id':['VDJDB0000000001'], 'gene':['TRB'], 'evidence_id':['motif_tcrnet:TRB'],
        'evidence_type':['motif_tcrnet'], 'evidence_source':['v1'], 'evidence_value':['cid'],
        'evidence_score':[5.], 'first_seen_release':['v1']}))
    monkeypatch.setattr(motifs, 'run_tcrnet', lambda *a,**k: ({}, ev))
    monkeypatch.setattr(motifs, 'per_epitope_report', lambda *a,**k: 0)
    for _ in range(2):
        motifs.run(tmp_path/'tables', tmp_path/'motifs', tmp_path/'reports', methods=('tcrnet',))
    evidence = vdjdb3.read_table(tmp_path/'tables','evidence')
    assert evidence.height == tables['evidence'].height + 1
    assert vdjdb3.read_table(tmp_path/'tables','vdjdb')['evidence.motif.tcrnet'].to_list() == [True,False]
    ev = ev.head(0)
    motifs.run(tmp_path/'tables', tmp_path/'motifs', tmp_path/'reports', methods=('tcrnet',))
    assert not vdjdb3.read_table(tmp_path/'tables','vdjdb')['evidence.motif.tcrnet'].any()


def test_chain_assembly_rejects_incomplete_partners():
    from vdjdb.assemble.master import CALL_FIELDS, FIX_FIELDS
    from vdjdb.assemble.tables import build_chains
    from vdjdb.schema import ALL_COLUMNS

    row = dict.fromkeys(ALL_COLUMNS, "")
    row.update(record_id="r1", species="HomoSapiens", TCR_hash="", **{
        "cdr3.alpha": "CAVRDSNYQLIW", "v.beta": "TRBV10-3*01",
        "d.beta": "TRBD1*01", "j.beta": "TRBJ2-7*01"})
    for gene in ("alpha", "beta"):
        for _, name, kind in FIX_FIELDS:
            row[f"{name}.{gene}"] = {str: "", bool: False, int: -1}[kind]
        for name, _ in CALL_FIELDS:
            row[f"{name}.{gene}"] = ""
        for segment in ("v", "j", "d"):
            row[f"__sub.{segment}.{gene}"] = ""
    with pytest.raises(ValueError, match="incomplete"):
        build_chains(pl.DataFrame([row]))


def test_explicit_coordinate_units_reach_both_export_formats(tmp_path, tables):
    values = {'v.end.nt': 4, 'v.end.aa': 2, 'j.start.nt': 8, 'j.start.aa': 2,
              'd.start.nt': 5, 'd.start.aa': 1, 'd.end.nt': 7, 'd.end.aa': 3}
    tables['chains'] = tables['chains'].with_columns(
        *[pl.lit(v, pl.Int64).alias(c) for c, v in values.items()])
    vdjdb3.write_all(tables, tmp_path)
    for name in ('chains', 'vdjdb'):
        disk = pl.read_parquet(tmp_path / f'{name}.parquet')
        for c, v in values.items():
            assert disk[c.replace('.', '_')].to_list() == [v, v]
    import io
    text = pl.read_csv(io.BytesIO(gzip.decompress((tmp_path / 'chains.tsv.gz').read_bytes())),
                       separator='\t')
    for c, v in values.items():
        assert text[c.replace('.', '_')].to_list() == [v, v]
    fields = {f['name']: f for f in json.loads((tmp_path / 'vdjdb.schema.json').read_text())['fields']}
    for c in values:
        assert fields[c]['dtype'] == 'Int64'
