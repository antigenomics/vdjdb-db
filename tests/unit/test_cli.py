"""The `vdjdb` command surface.

`cli.py` had no test at all: 173 statements, 0 % covered, while being the only way anyone runs this
build. These cover what can be exercised without a full build - argument handling, the projections
of the field registry, and the exit codes callers and workflows depend on.

The expensive subcommands (`build`, `motifs`, `summary`, `refs`) are represented by their argument
validation only. `build.yml` runs the real thing.
"""
from __future__ import annotations

import json

import pytest
from typer.testing import CliRunner

from vdjdb import __version__
from vdjdb.cli import app
from vdjdb.schema import TABLES, header, render_meta, render_slim_meta

runner = CliRunner()

#: Every command `vdjdb --help` must offer. A command disappearing from the app is a packaging
#: break that no other test would notice.
COMMANDS = ["version", "qc", "schema", "build", "make", "rules", "convert", "motifs", "summary",
            "refs", "release", "changelog", "diff"]


def test_help_lists_every_command():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    missing = [c for c in COMMANDS if c not in result.output]
    assert not missing, f"missing from --help: {missing}"


@pytest.mark.parametrize("command", COMMANDS)
def test_each_command_has_help(command):
    """A command whose options fail to construct raises here rather than in a workflow."""
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0, result.output


def test_version_reports_package_and_root():
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert __version__ in result.output
    assert "root " in result.output


@pytest.mark.parametrize("table", sorted(TABLES))
def test_schema_header_matches_registry(table):
    """The printed header is the registry's, for every declared table."""
    result = runner.invoke(app, ["schema", "--table", table, "--format", "header"])
    assert result.exit_code == 0, result.output
    assert result.output.rstrip("\n") == header(table)


def test_schema_meta_matches_renderer():
    result = runner.invoke(app, ["schema", "--table", "vdjdb", "--format", "meta"])
    assert result.exit_code == 0
    assert result.output == render_meta("vdjdb")


def test_schema_meta_uses_the_slim_renderer_for_slim():
    """slim's metadata is two columns; every other table's is the eight-column form."""
    result = runner.invoke(app, ["schema", "--table", "slim", "--format", "meta"])
    assert result.exit_code == 0
    assert result.output == render_slim_meta("slim")
    assert result.output != render_meta("slim")


def test_schema_json_is_parseable_and_one_entry_per_column():
    result = runner.invoke(app, ["schema", "--table", "records", "--format", "json"])
    assert result.exit_code == 0
    payload = json.loads(result.output)
    assert isinstance(payload, list) and payload
    assert len(payload) == len(header("records").split("\t"))


def test_schema_rejects_an_unknown_table():
    result = runner.invoke(app, ["schema", "--table", "no-such-table"])
    assert result.exit_code == 2
    assert "unknown table" in result.output


def test_schema_rejects_an_unknown_format():
    result = runner.invoke(app, ["schema", "--table", "vdjdb", "--format", "yaml"])
    assert result.exit_code == 2
    assert "unknown format" in result.output


def test_make_rejects_an_unknown_projection():
    result = runner.invoke(app, ["make", "vdjdb4"])
    assert result.exit_code == 2
    assert "unknown projection" in result.output


def test_convert_rejects_an_unknown_target():
    result = runner.invoke(app, ["convert", "parquet"])
    assert result.exit_code == 2
    assert "unknown target" in result.output


@pytest.mark.parametrize("args", [
    [],                                     # neither source
    ["--tables", "out/tables", "--legacy", "out/legacy/vdjdb.txt"],   # both
])
def test_convert_airr_needs_exactly_one_source(args, tmp_path):
    result = runner.invoke(app, ["convert", "airr", "--out", str(tmp_path), *args])
    assert result.exit_code == 2
    assert "exactly one" in result.output


def test_qc_accepts_a_clean_chunk(tmp_path):
    """A well-formed chunk exits 0, which is what a submission pull request depends on."""
    chunk = tmp_path / "PMID_1.txt"
    chunk.write_text(
        "cdr3.alpha\tv.alpha\tj.alpha\tcdr3.beta\tv.beta\td.beta\tj.beta\tspecies\tmhc.a\tmhc.b\t"
        "mhc.class\tantigen.epitope\tantigen.gene\tantigen.species\treference.id\t"
        "method.identification\tmethod.frequency\tmethod.singlecell\tmethod.sequencing\t"
        "method.verification\tmeta.study.id\tmeta.cell.subset\tmeta.subject.cohort\t"
        "meta.subject.id\tmeta.replica.id\tmeta.clone.id\tmeta.epitope.id\tmeta.tissue\t"
        "meta.donor.MHC\tmeta.donor.MHC.method\tmeta.structure.id\n")
    result = runner.invoke(app, ["qc", str(chunk)])
    assert result.exit_code == 0, result.output


def test_qc_strict_fails_on_an_unreadable_chunk(tmp_path):
    bad = tmp_path / "PMID_2.txt"
    bad.write_text("not\ta\tchunk\n1\t2\t3\n")
    result = runner.invoke(app, ["qc", str(bad)])
    assert result.exit_code != 0


def test_qc_writes_a_report(tmp_path):
    bad = tmp_path / "PMID_3.txt"
    bad.write_text("not\ta\tchunk\n1\t2\t3\n")
    report = tmp_path / "qc.tsv"
    runner.invoke(app, ["qc", str(bad), "--report", str(report)])
    assert report.exists()
    assert report.read_text().startswith("file\trow\tcode\tdetail")


@pytest.mark.parametrize("tag", ["2026.09.1", "v2026.9.1", "v2026-09-01", "v2026.09", "latest"])
def test_release_rejects_a_tag_that_is_not_the_declared_scheme(tag, tmp_path):
    result = runner.invoke(app, ["release", "--tag", tag, "--build", str(tmp_path)])
    assert result.exit_code != 0


def test_diff_reports_a_missing_reference(tmp_path):
    result = runner.invoke(app, ["diff", str(tmp_path / "absent.zip"), str(tmp_path)])
    assert result.exit_code != 0
