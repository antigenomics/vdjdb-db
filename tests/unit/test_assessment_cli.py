"""Selected assessments recompute one stage without invoking the full assembly line."""
from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from vdjdb.cli import app
from vdjdb.emit.vdjdb3 import read_table

runner = CliRunner()


def observations():
    return pl.DataFrame({
        "antigen.epitope": ["AAAAAAAAA"], "antigen.species": ["HomoSapiens"],
        "antigen.gene": [""], "species": ["HomoSapiens"], "mhc.class": ["MHCI"],
        "mhc.a": ["HLA-A*02:01"], "mhc.b": ["B2M"], "reference.id": ["PMID:1"],
    })


def test_selected_chunks_write_only_assessment_and_pass_resource_options(tmp_path, monkeypatch):
    from vdjdb.assemble import assessment, master
    calls = []
    monkeypatch.setattr(master, "build_master", lambda paths: calls.append(paths) or observations())
    real = assessment.build_assessment

    def assess(records, **kwargs):
        assert kwargs == {"reference": None, "jobs": 2}
        return real(records, **kwargs)

    monkeypatch.setattr(assessment, "build_assessment", assess)
    result = runner.invoke(app, ['assess-epitopes', 'chunks/PMID_1.tsv', '--jobs', '2',
                                 '--out', str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert calls == [[Path('chunks/PMID_1.tsv')]]
    assert {p.name for p in tmp_path.iterdir()} == {
        'epitope_assessment.parquet', 'epitope_assessment.tsv', 'assessment-timings.tsv'}
    got = read_table(tmp_path, 'epitope_assessment')
    assert got['assessment.status'].to_list() == ['reference_not_supplied']
    assert got['records'].to_list() == [1]


def test_table_source_recomputes_from_records_without_reading_old_assessment(tmp_path):
    from vdjdb.schema import TIDY_NAMES
    observations().rename(TIDY_NAMES, strict=False).write_parquet(tmp_path / 'records.parquet')
    (tmp_path / 'epitope_assessment.parquet').write_text('invalid stale result')
    out = tmp_path / 'new'
    result = runner.invoke(app, ['assess-epitopes', '--tables', str(tmp_path), '--out', str(out)])
    assert result.exit_code == 0, result.output
    assert read_table(out, 'epitope_assessment')['records'].sum() == 1


def test_assessment_requires_one_explicit_source(monkeypatch):
    from vdjdb.assemble import master
    from vdjdb.emit import vdjdb3

    def unexpected_read(*args, **kwargs):
        pytest.fail('ambiguous input must be rejected before reading any source')

    monkeypatch.setattr(master, 'build_master', unexpected_read)
    monkeypatch.setattr(vdjdb3, 'read_table', unexpected_read)
    for args in [[], ['chunks/PMID_1.tsv', '--tables', 'out/tables']]:
        result = runner.invoke(app, ['assess-epitopes', *args])
        assert result.exit_code == 2
