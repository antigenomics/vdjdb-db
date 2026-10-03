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


#: A row that passes every rule. Each of the three tests below used to pass a header with no data
#: rows, which `no-data-rows` now reports: `ChunkQC.check_exist` raised `ValueError("Empty file")` on
#: exactly that shape, and a chunk contributing no records was reaching every gate green.
CLEAN = {"cdr3.beta": "CASSIRSSYEQYF", "v.beta": "TRBV10-3*01", "j.beta": "TRBJ2-7*01",
         "species": "HomoSapiens", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M", "mhc.class": "MHCI",
         "antigen.epitope": "GILGFVFTL", "antigen.gene": "M", "antigen.species": "InfluenzaA",
         "reference.id": "PMID:1"}


def clean_row(**over):
    return blank_row(**(CLEAN | over))


def test_a_clean_chunk_passes(tmp_path):
    assert run_qc([chunk(tmp_path, rows=[clean_row()])], strict=True) == 0


def test_a_chunk_with_no_records_is_fatal(tmp_path):
    """`check_exist` raised on this and nothing here reported it: `empty` needs a file with no lines
    at all, and the reader returns a 0-row frame without complaint."""
    assert "no-data-rows" not in ADVISORY
    assert run_qc([chunk(tmp_path)], strict=True) == 1


def test_crlf_is_reported_but_never_fatal(tmp_path):
    """Advisory by design: a hard gate would have blocked every unrelated submission until the
    line-ending normalisation landed."""
    p = tmp_path / "PMID_2.txt"
    p.write_bytes((HEADER + "\r\n" + clean_row() + "\r\n").encode())
    assert "crlf" in ADVISORY
    assert run_qc([p], strict=True) == 0


def test_an_unknown_column_is_advisory(tmp_path):
    p = chunk(tmp_path, rows=[clean_row() + "\tsomething"], extra_header="\tnot.a.column")
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


# --- a named segment whose chain has no CDR3 ----------------------------------------------------

def test_a_segment_call_with_no_cdr3_is_reported_and_is_advisory():
    """The call is information, the row is kept, but the chain cannot reach any output.

    Every shipped table is keyed on the CDR3, so a V or J named without one is carried in `chunks/`
    and dropped by the build. A submitter should hear that while they can still supply the sequence.
    """
    from vdjdb.qc.rules import RULES

    rows = pl.DataFrame({
        "chunk.file": ["c.txt"] * 4,
        "chunk.row": [0, 1, 2, 3],
        # a V named with no alpha CDR3; a J likewise; a complete chain; a chain absent entirely
        "cdr3.alpha": ["", "", "CAVRDSNYQLIW", ""],
        "v.alpha": ["TRAV12-2*01", "", "TRAV12-2*01", ""],
        "j.alpha": ["", "TRAJ33*01", "TRAJ33*01", ""],
        "cdr3.beta": ["CASSIRSSYEQYF"] * 4,
        "v.beta": ["TRBV10-3*01"] * 4,
        "j.beta": ["TRBJ2-7*01"] * 4,
    })
    # one rule, not `check`, so the fixture does not have to carry every column every rule reads
    failed = rows.filter(~RULES["segment call with no cdr3"])
    assert sorted(failed["chunk.row"].to_list()) == [0, 1], "only the two incomplete alpha chains"
    assert "segment call with no cdr3" in ADVISORY, "a chain that cannot ship is not a broken row"


# --- a structure id that is not a PDB entry id ---------------------------------------------------

def test_a_structure_id_that_is_not_a_pdb_id_is_reported_and_is_advisory():
    """The field awards the top confidence score, so anything but a PDB id awards it for nothing.

    `score.confidence` gives 3 outright when `meta.structure.id` is non-empty, above every
    sequencing and specificity term, because a solved TCR:pMHC complex is direct proof of binding.
    2,765 rows in the corpus hold a figure or table reference there instead.
    """
    from vdjdb.qc.rules import RULES

    rows = pl.DataFrame({
        "chunk.file": ["c.txt"] * 6,
        "chunk.row": [0, 1, 2, 3, 4, 5],
        "meta.structure.id": [
            "1AO7",                               # a PDB id
            "6uon",                               # lower case is still a PDB id
            "",                                   # blank is the normal case
            "Fig 9, Supp Fig 5, Supp Table 5-8",  # the 400-row case in menon_etal_2024.tsv
            "56I",                                # three characters, the PMID_28423320.tsv case
            "ABCD",                               # four alphanumerics, but a PDB id starts with a digit
        ],
    })
    failed = rows.filter(~RULES["structure id is not a PDB id"])
    assert sorted(failed["chunk.row"].to_list()) == [3, 4, 5]
    assert "structure id is not a PDB id" in ADVISORY, (
        "2,765 corpus rows fail it; what to do with them is a curation decision")


@pytest.mark.parametrize("values,rule", [
    ({"cdr3.alpha": "", "cdr3.beta": ""}, "no.cdr3"),
    ({"antigen.epitope": ""}, "no.antigen.seq"),
    ({"mhc.a": ""}, "no.mhc"),
    ({"mhc.b": ""}, "no.mhc"),
    ({"mhc.class": ""}, "bad mhc.class"),
    ({"mhc.b": "HLA-DRB1*01:01"}, "mhc class/partner mismatch"),
    ({"mhc.class": "MHCII"}, "mhc class/partner mismatch"),
    ({"mhc.a": "B2M"}, "mhc class/partner mismatch"),
])
def test_incomplete_observations_fail_qc_and_direct_build(tmp_path, values, rule):
    from vdjdb.assemble.master import build_master

    path = chunk(tmp_path, rows=[clean_row(**values)])
    assert rule not in ADVISORY
    assert run_qc([path], strict=True) == 1
    with pytest.raises(ValueError, match=rule):
        build_master([path])


@pytest.mark.parametrize("chain", ["alpha", "beta"])
def test_one_junction_and_missing_segment_calls_are_allowed(tmp_path, chain):
    from vdjdb.io.chunks import read_chunk
    from vdjdb.qc.rules import assert_complete

    values = {"cdr3.beta": "", "v.beta": "", "j.beta": "",
              f"cdr3.{chain}": "CASSIRSSYEQYF"}
    path = chunk(tmp_path, rows=[clean_row(**values)])
    assert run_qc([path], strict=True) == 0
    assert_complete(read_chunk(path))
