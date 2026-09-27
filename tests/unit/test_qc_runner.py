"""`vdjdb qc`: what fails a submission, what is only reported, and what the report contains.

`qc/runner.py` was 28 % covered while being the gate every chunk pull request hits. The split
between advisory and fatal findings is the part that matters: 99 of 230 chunks are CRLF, so a hard
gate on that would block every unrelated submission, while a bad epitope has to fail.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.qc.lint import Finding
from vdjdb.qc.runner import ADVISORY, _report_frame, run_qc
from vdjdb.schema import ALL_COLUMNS

#: A chunk header the reader accepts, built from the declared column set so it cannot drift.
HEADER = "\t".join(ALL_COLUMNS)


def chunk(tmp_path, name="PMID_1.txt", rows=(), extra_header=""):
    p = tmp_path / name
    head = HEADER + extra_header
    p.write_text(head + "\n" + "".join(r + "\n" for r in rows))
    return p


def blank_row(**over):
    cells = {c: "" for c in ALL_COLUMNS}
    cells.update(over)
    return "\t".join(cells[c] for c in ALL_COLUMNS)


def test_a_clean_chunk_passes(tmp_path):
    assert run_qc([chunk(tmp_path)], strict=True) == 0


def test_crlf_is_reported_but_never_fatal(tmp_path):
    """Advisory by design: 99 of 230 chunks are CRLF until the .tsv migration lands."""
    p = tmp_path / "PMID_2.txt"
    p.write_bytes((HEADER + "\r\n").encode())
    assert "crlf" in ADVISORY
    assert run_qc([p], strict=True) == 0


def test_an_unknown_column_is_advisory(tmp_path):
    p = chunk(tmp_path, extra_header="\tnot.a.column")
    assert "unknown-column" in ADVISORY
    assert run_qc([p], strict=True) == 0


def test_strict_fails_and_no_strict_does_not(tmp_path):
    """The same findings, two exit codes. `--no-strict` is how a report is produced anyway."""
    p = chunk(tmp_path, rows=[blank_row(**{"antigen.epitope": "not an epitope!",
                                           "species": "HomoSapiens"})])
    assert run_qc([p], strict=True) == 1
    assert run_qc([p], strict=False) == 0


def test_a_report_is_written_even_when_the_chunk_cannot_be_parsed(tmp_path):
    """The header being wrong is the commonest submission mistake, and it used to raise.

    chunk-check uploads this file and its pull-request comment is built from it, so a run that
    produces no report reports "no findings" for the case that needed one most.
    """
    bad = tmp_path / "PMID_3.txt"
    bad.write_text("not\ta\tchunk\n1\t2\t3\n")
    report = tmp_path / "qc.tsv"
    assert run_qc([bad], strict=True, report=report) == 1
    body = report.read_text()
    assert body.startswith("file\trow\tcode\tdetail")
    assert "unreadable" in body


def test_the_report_carries_both_tiers(tmp_path):
    p = chunk(tmp_path, extra_header="\tnot.a.column",
              rows=[blank_row(**{"antigen.epitope": "not an epitope!", "species": "HomoSapiens",
                                 "not.a.column": ""})])
    report = tmp_path / "qc.tsv"
    run_qc([p], strict=False, report=report)
    codes = set(pl.read_csv(report, separator="\t")["code"])
    assert "unknown-column" in codes, "a text-level finding"
    assert len(codes) > 1, codes


def test_the_report_is_sorted_so_two_runs_agree(tmp_path):
    p = chunk(tmp_path, extra_header="\tzzz\taaa")
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    run_qc([p], strict=False, report=a)
    run_qc([p], strict=False, report=b)
    assert a.read_text() == b.read_text()


def test_report_frame_has_the_declared_columns_when_empty():
    out = _report_frame([], pl.DataFrame(schema={"chunk.file": pl.Utf8, "chunk.row": pl.UInt32,
                                                 "rule": pl.Utf8}))
    assert out.columns == ["file", "row", "code", "detail"]
    assert out.height == 0


def test_report_frame_keeps_the_detail_of_a_text_finding():
    out = _report_frame([Finding(file="PMID_9.txt", code="crlf", detail="line 3")],
                        pl.DataFrame(schema={"chunk.file": pl.Utf8, "chunk.row": pl.UInt32,
                                             "rule": pl.Utf8}))
    assert out.row(0, named=True) == {"file": "PMID_9.txt", "row": 0, "code": "crlf",
                                      "detail": "line 3"}


@pytest.mark.parametrize("code", sorted(ADVISORY))
def test_every_advisory_code_is_one_something_can_emit(code):
    """A name left in ADVISORY after its rule was renamed downgrades nothing, and says nothing."""
    import re

    from vdjdb.qc.rules import RULES

    # `check` emits every RULES key plus `duplicate`, which it computes inline as `__duplicate`.
    emitted = set(RULES) | {"duplicate"}
    # The linter's codes are the second argument of each `Finding(...)` it constructs.
    emitted |= set(re.findall(r'Finding\([^,]+,\s*"([^"]+)"',
                              Path("src/vdjdb/qc/lint.py").read_text()))
    assert code in emitted, f"{code} is in ADVISORY but nothing emits it; emitted: {sorted(emitted)}"
