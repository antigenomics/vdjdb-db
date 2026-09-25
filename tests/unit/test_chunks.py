"""The chunk reader and the vectorised row rules."""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.io.chunks import PROVENANCE, READABLE, dedup, independently_reported, read_chunk, read_chunks
from vdjdb.qc.rules import RULES, check
from vdjdb.schema import ALL_COLUMNS, CHUNK_DEDUP_KEY


def _chunk(tmp: Path, name: str, rows: list[dict[str, str]], *,
           extra: tuple[str, ...] = (), sep: str = "\n", prefix: str = "") -> Path:
    cols = list(ALL_COLUMNS) + list(extra)
    lines = ["\t".join(cols)]
    lines += ["\t".join(r.get(c, "") for c in cols) for r in rows]
    p = tmp / name
    p.write_text(prefix + sep.join(lines) + sep, encoding="utf-8")
    return p


def _row(**kw: str) -> dict[str, str]:
    base = {"cdr3.beta": "CASSA", "v.beta": "TRBV1*01", "j.beta": "TRBJ1*01",
            "species": "HomoSapiens", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M",
            "mhc.class": "MHCI", "antigen.epitope": "GILGFVFTL", "antigen.gene": "M",
            "antigen.species": "InfluenzaA", "reference.id": "PMID:1"}
    return base | kw


# --------------------------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------------------------

def test_every_cell_is_a_string_and_every_gap_is_an_empty_string(tmp_path: Path) -> None:
    """One missing marker. The pandas None/NaN/"" ambiguity shipped real bugs."""
    df = read_chunk(_chunk(tmp_path, "a.txt", [_row()]))
    assert df.null_count().sum_horizontal().item() == 0
    assert set(df.schema.values()) <= {pl.Utf8, pl.UInt32}
    assert df["cdr3.alpha"][0] == ""


def test_provenance_is_carried(tmp_path: Path) -> None:
    df = read_chunk(_chunk(tmp_path, "a.txt", [_row(), _row(**{"meta.subject.id": "d2"})]))
    assert df["chunk.file"].to_list() == ["a.txt", "a.txt"]
    assert df["chunk.row"].to_list() == [1, 2], "1-based, counting data rows"


def test_columns_are_the_registry_order_plus_provenance(tmp_path: Path) -> None:
    df = read_chunk(_chunk(tmp_path, "a.txt", [_row()]))
    assert df.columns == [*READABLE, *PROVENANCE]


def test_crlf_and_bom_are_tolerated(tmp_path: Path) -> None:
    """99 of 230 chunks were CRLF; a BOM appears on spreadsheet exports."""
    p = _chunk(tmp_path, "a.txt", [_row()], sep="\r\n", prefix="﻿")
    df = read_chunk(p)
    assert df.height == 1
    assert df["reference.id"][0] == "PMID:1", "no stray carriage return on the last column"


def test_the_one_capitalised_column_name_is_accepted(tmp_path: Path) -> None:
    p = _chunk(tmp_path, "a.txt", [_row()], extra=("Comment",))
    assert read_chunk(p)["comment"][0] == ""


def test_an_unknown_column_is_dropped_not_carried(tmp_path: Path) -> None:
    p = _chunk(tmp_path, "a.txt", [_row()], extra=("nonsense",))
    assert "nonsense" not in read_chunk(p).columns


def test_a_missing_required_column_is_refused(tmp_path: Path) -> None:
    cols = [c for c in ALL_COLUMNS if c != "mhc.class"]
    p = tmp_path / "a.txt"
    p.write_text("\t".join(cols) + "\n" + "\t".join("" for _ in cols) + "\n")
    with pytest.raises(ValueError, match="missing required columns"):
        read_chunk(p)


def test_a_duplicate_column_is_refused_rather_than_guessed(tmp_path: Path) -> None:
    cols = [*ALL_COLUMNS, "species"]
    p = tmp_path / "a.txt"
    p.write_text("\t".join(cols) + "\n" + "\t".join("" for _ in cols) + "\n")
    with pytest.raises(ValueError, match="duplicate columns"):
        read_chunk(p)


def test_files_are_read_in_sorted_order(tmp_path: Path) -> None:
    """``os.listdir`` order is why the release cannot be reproduced byte-for-byte."""
    for n in ("c.txt", "a.txt", "b.txt"):
        _chunk(tmp_path, n, [_row(**{"meta.study.id": n})])
    df = read_chunks(sorted(tmp_path.glob("*.txt")), deduplicate=False)
    assert df["chunk.file"].to_list() == ["a.txt", "b.txt", "c.txt"]


def test_reading_is_reproducible(tmp_path: Path) -> None:
    files = [_chunk(tmp_path, f"{n}.txt", [_row(**{"meta.subject.id": f"d{i}"})
                                           for i in range(3)]) for n in "abc"]
    digests = {read_chunks(files).hash_rows().sum() for _ in range(5)}
    assert len(digests) == 1


# --------------------------------------------------------------------------------------------
# Deduplication
# --------------------------------------------------------------------------------------------

def test_dedup_is_per_chunk_because_two_chunks_are_two_papers(tmp_path: Path) -> None:
    """Within a chunk, a repeated row is a duplicate. Across chunks it is independent replication.

    Collapsing the cross-chunk pair would delete the strongest evidence the database carries.
    """
    a = _chunk(tmp_path, "a.txt", [_row(), _row()])
    b = _chunk(tmp_path, "b.txt", [_row()])
    df = read_chunks([a, b])
    assert df.height == 2, "the within-chunk duplicate goes, the second paper's report stays"
    assert independently_reported(df).height == 1


def test_dedup_keeps_the_first_occurrence_and_the_order(tmp_path: Path) -> None:
    p = _chunk(tmp_path, "a.txt", [_row(**{"meta.tissue": "PBMC"}),
                                   _row(**{"meta.tissue": "PBMC"}),
                                   _row(**{"meta.tissue": "LN"})])
    df = dedup(read_chunk(p))
    assert df["meta.tissue"].to_list() == ["PBMC", "LN"]
    assert df["chunk.row"].to_list() == [1, 3]


def test_a_field_outside_the_dedup_key_does_not_split_a_duplicate(tmp_path: Path) -> None:
    """``meta.structure.id`` is not identity; two rows differing only there are one record."""
    assert "meta.structure.id" not in CHUNK_DEDUP_KEY
    p = _chunk(tmp_path, "a.txt", [_row(**{"meta.structure.id": "5EUO"}),
                                   _row(**{"meta.structure.id": "5euo"})])
    assert dedup(read_chunk(p)).height == 1


# --------------------------------------------------------------------------------------------
# Row rules
# --------------------------------------------------------------------------------------------

def _findings(tmp: Path, row: dict[str, str]) -> set[str]:
    df = read_chunk(_chunk(tmp, "a.txt", [row]))
    return set(check(df)["rule"].to_list())


def test_a_clean_row_raises_nothing(tmp_path: Path) -> None:
    assert _findings(tmp_path, _row()) == set()


@pytest.mark.parametrize(("field", "value", "rule"), [
    ("cdr3.beta", "CASSX", "bad cdr3.beta"),
    ("cdr3.beta", "CAS", "bad cdr3.beta"),
    ("v.beta", "TRAV1*01", "bad v.beta"),
    ("j.beta", "TRBV1*01", "bad j.beta"),
    ("species", "Human", "bad species"),
    ("mhc.class", "MHC1", "bad mhc.class"),
    ("mhc.a", "HLA-A2", "bad mhc.a"),
    ("antigen.epitope", "GILGFVFTLX", "bad antigen.epitope"),
    ("reference.id", "see figure 3", "bad reference.id"),
])
def test_each_rule_fires_on_its_own_defect(tmp_path: Path, field: str, value: str,
                                           rule: str) -> None:
    assert rule in _findings(tmp_path, _row(**{field: value}))


def test_a_row_with_no_cdr3_at_all_is_caught(tmp_path: Path) -> None:
    assert "no.cdr3" in _findings(tmp_path, _row(**{"cdr3.beta": ""}))


def test_a_row_with_no_epitope_is_caught(tmp_path: Path) -> None:
    assert "no.antigen.seq" in _findings(tmp_path, _row(**{"antigen.epitope": ""}))


def test_a_row_with_half_an_mhc_is_caught(tmp_path: Path) -> None:
    assert "no.mhc" in _findings(tmp_path, _row(**{"mhc.b": ""}))


def test_murine_mhc_is_deliberately_unchecked(tmp_path: Path) -> None:
    """``is_MHC_valid`` only constrains HLA spellings, which is why the murine MHC-II
    fragmentation (``I-Ab`` 3,274 vs ``H2-IAb`` 113 vs ``H2-Ab1`` 9) never tripped QC.

    Changing it here would fail 230 chunks at once; phase 9's rule tables are where it is fixed.
    """
    assert "bad mhc.a" not in _findings(tmp_path, _row(**{"mhc.a": "I-Ab", "mhc.b": "H2-Ab1",
                                                          "species": "MusMusculus"}))


def test_every_rule_is_expressed_as_an_invariant(tmp_path: Path) -> None:
    """Rules read as what they protect, not as what they catch, so a clean row satisfies all."""
    df = read_chunk(_chunk(tmp_path, "a.txt", [_row()]))
    for name, expr in RULES.items():
        assert df.select(expr).to_series().all(), f"{name} rejects a clean row"


def test_duplicate_is_reported_per_chunk(tmp_path: Path) -> None:
    df = read_chunk(_chunk(tmp_path, "a.txt", [_row(), _row(), _row()]))
    dupes = check(df).filter(pl.col("rule") == "duplicate")
    assert dupes["chunk.row"].to_list() == [2, 3], "the first occurrence is not a duplicate"
