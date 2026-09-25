"""Chunk lint checks. These are the rules a submitter hits in CI, so they get real coverage."""
from __future__ import annotations

import pytest

from vdjdb.qc.lint import lint_file
from vdjdb.schema.tables import ALL_COLUMNS, CHUNK_DEDUP_KEY, COMPLEX_COLUMNS, TOLERATED_DROPPED

HEADER = "\t".join(ALL_COLUMNS)
ROW = "\t".join(
    "CASSL" if c.startswith("cdr3") else "HomoSapiens" if c == "species"
    else "PMID:12345678" if c == "reference.id" else ""
    for c in ALL_COLUMNS
)


def write(tmp_path, text: str, *, name: str = "PMID_1.txt", newline: str = "\n"):
    p = tmp_path / name
    p.write_bytes(text.replace("\n", newline).encode("utf-8"))
    return p


def codes(findings) -> set[str]:
    return {f.code for f in findings}


def test_clean_chunk_has_no_findings(tmp_path):
    assert lint_file(write(tmp_path, f"{HEADER}\n{ROW}\n")) == []


def test_crlf_is_flagged(tmp_path):
    assert "crlf" in codes(lint_file(write(tmp_path, f"{HEADER}\n{ROW}\n", newline="\r\n")))


def test_bom_is_flagged(tmp_path):
    p = tmp_path / "PMID_2.txt"
    p.write_bytes(b"\xef\xbb\xbf" + f"{HEADER}\n{ROW}\n".encode())
    assert "bom" in codes(lint_file(p))


def test_empty_column_name_is_flagged(tmp_path):
    # The real failure mode: a pandas index artefact leaves a bare leading tab.
    assert "empty-column-name" in codes(lint_file(write(tmp_path, f"\t{HEADER}\n\t{ROW}\n")))


def test_duplicate_column_name_is_flagged(tmp_path):
    assert "duplicate-column-name" in codes(lint_file(write(tmp_path, f"{HEADER}\tspecies\n{ROW}\t\n")))


def test_prose_used_as_column_name_is_flagged(tmp_path):
    prose = "optional columns that ensure correct confidence ranking of a given entry"
    assert "prose-column-name" in codes(lint_file(write(tmp_path, f"{HEADER}\t{prose}\n{ROW}\t\n")))


def test_missing_required_column_is_flagged(tmp_path):
    reduced = "\t".join(c for c in ALL_COLUMNS if c != "mhc.class")
    row = "\t".join("" for _ in range(len(ALL_COLUMNS) - 1))
    found = [f for f in lint_file(write(tmp_path, f"{reduced}\n{row}\n")) if f.code == "missing-required-column"]
    assert found and "mhc.class" in found[0].detail


@pytest.mark.parametrize("extra", sorted(TOLERATED_DROPPED))
def test_tolerated_columns_are_not_unknown(tmp_path, extra):
    """The five columns the legacy build silently discards must not be reported as unknown.

    Listing them is what makes the drop deliberate rather than accidental -- see ROADMAP section 9.
    """
    found = lint_file(write(tmp_path, f"{HEADER}\t{extra}\n{ROW}\t\n"))
    assert "unknown-column" not in codes(found)


def test_genuinely_unknown_column_is_flagged(tmp_path):
    assert "unknown-column" in codes(lint_file(write(tmp_path, f"{HEADER}\tnot.a.column\n{ROW}\t\n")))


def test_bad_reference_id_form_is_flagged(tmp_path):
    bad = ROW.replace("PMID:12345678", "Smith et al 2020")
    assert "reference-id-form" in codes(lint_file(write(tmp_path, f"{HEADER}\n{bad}\n")))


def test_schema_tuples_are_self_consistent():
    assert len(ALL_COLUMNS) == len(set(ALL_COLUMNS)), "duplicate column in ALL_COLUMNS"
    assert set(COMPLEX_COLUMNS) <= set(ALL_COLUMNS)
    assert set(CHUNK_DEDUP_KEY) <= set(ALL_COLUMNS)
    assert not (TOLERATED_DROPPED & set(ALL_COLUMNS)), "a tolerated-dropped column is also a real one"
