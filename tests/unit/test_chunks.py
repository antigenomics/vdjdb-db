"""The chunk reader and the vectorised row rules."""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.config import Paths
from vdjdb.io.chunks import (
    PROVENANCE,
    READABLE,
    chunk_files,
    dedup,
    independently_reported,
    read_chunk,
    read_chunks,
)
from vdjdb.qc.rules import RULES, check
from vdjdb.schema import ALL_COLUMNS


def _chunk(tmp: Path, name: str, rows: list[dict[str, str]], *,
           extra: tuple[str, ...] = (), sep: str = "\n", prefix: str = "") -> Path:
    cols = list(ALL_COLUMNS) + list(extra)
    lines = ["\t".join(cols)]
    lines += ["\t".join(r.get(c, "") for c in cols) for r in rows]
    p = tmp / name
    p.write_text(prefix + sep.join(lines) + sep, encoding="utf-8")
    return p


def _row(**kw: str) -> dict[str, str]:
    # Real, functional human calls. `TRBV1*01` stood here until #634 made IMGT's F / ORF / P verdict
    # a rule: it is a human *pseudogene*, so the "clean" row was never clean and
    # `test_every_rule_is_expressed_as_an_invariant` failed on a rule that was working.
    base = {"cdr3.beta": "CASSA", "v.beta": "TRBV10-3*01", "j.beta": "TRBJ2-7*01",
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


@pytest.mark.parametrize("field", ["method.identification", "method.verification",
    "method.frequency", "meta.structure.id", "meta.donor.MHC", "meta.epitope.id"])
def test_observation_metadata_distinguishes_rows(tmp_path: Path, field: str) -> None:
    p = _chunk(tmp_path, "a.txt", [_row(**{field: "one"}), _row(**{field: "two"})])
    assert dedup(read_chunk(p)).height == 2


def test_the_read_path_loses_nothing_but_declared_duplicates() -> None:
    """Every row `chunks/` holds is a record or a within-chunk duplicate on `CHUNK_DEDUP_KEY`.

    `EXPECTED_RECORDS` in `tests/release/test_tables_contract.py` pins the total, and a total tells
    you a number moved without telling you why. This pins the *reason*: the only transformation
    between the files and the frame the build reads is the declared deduplication. A filter added to
    the read path - dropping a blank CDR3, skipping a species with no germline - would lose curated
    records silently, and the release tier is where it would otherwise be noticed, needing a
    reference zip a laptop does not have.

    Runs on the corpus rather than a fixture, because a fixture cannot catch a filter conditioned on
    a value only real data carries.
    """
    # Counted off disk, because that is the only count no part of the read path can influence.
    # Comparing two reads against each other cannot catch a filter inside `read_chunk`: it would
    # apply to both sides of the comparison and the equality would still hold.
    on_disk = sum(len([ln for ln in f.read_bytes().split(b"\n") if ln.strip()]) - 1
                  for f in chunk_files())
    raw = read_chunks(deduplicate=False)
    assert raw.height == on_disk, (
        f"{on_disk - raw.height} curated line(s) did not survive reading. The reader normalises "
        "and annotates; it must never select.")

    kept = read_chunks()
    assert kept.equals(dedup(raw)), "the read path applied something other than dedup()"
    # Provenance survives for every kept row, so any record traces back to a curator's line.
    assert kept.select("chunk.file", "chunk.row").n_unique() == kept.height


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


#: Input directories holding chunks the build must never read. `pending/` is the current format
#: waiting on a reference the build lacks; `withheld/` predates the specification. See `CLAUDE.md`.
QUARANTINED = ("pending", "withheld")


@pytest.mark.parametrize("directory", QUARANTINED)
def test_a_quarantined_chunk_cannot_reach_the_build(directory: str) -> None:
    """The only thing keeping these out is that `chunk_files` reads one directory, not a tree.

    Nothing else in the repository names them, so a later `rglob("PMID_*.txt")` would pull them in
    silently -- and `pending/PMID_22058411.txt` would enter the build carrying a species no part of
    it can handle. This test is what fails if that happens.
    """
    root = Paths.discover().root
    quarantined = root / directory
    assert quarantined.is_dir(), f"{directory}/ is missing"
    assert list(quarantined.glob("*.txt")), f"{directory}/ holds no chunks, so this proves nothing"

    read = {p.resolve() for p in chunk_files()}
    assert read, "no chunks were read at all"
    intruders = sorted(p.name for p in quarantined.glob("*.txt") if p.resolve() in read)
    assert not intruders, f"the build reads {len(intruders)} file(s) from {directory}/: {intruders}"


def test_the_two_quarantine_directories_hold_different_formats() -> None:
    """The rule that decides which directory a blocked chunk goes to, asserted on the real files.

    `pending/` is for a chunk whose header matches a shipping chunk; `withheld/` is for one that
    predates it. Get this backwards and a curator is told to re-export a file that needs no export.
    """
    root = Paths.discover().root
    shipping = len(chunk_files()[0].read_text(encoding="utf-8").split("\n", 1)[0].split("\t"))

    for p in sorted((root / "pending").glob("*.txt")):
        got = len(p.read_text(encoding="utf-8").split("\n", 1)[0].split("\t"))
        assert got == shipping, (
            f"pending/{p.name} has {got} columns against the current {shipping}. A chunk the reader "
            "cannot parse belongs in withheld/, which says the file must be re-exported.")

    for p in sorted((root / "withheld").glob("*.txt")):
        got = len(p.read_text(encoding="utf-8").split("\n", 1)[0].split("\t"))
        assert got != shipping, (
            f"withheld/{p.name} already has the current {shipping}-column header, so it parses. It "
            "belongs in pending/, which says the build is what has to change.")


def test_compressed_chunk_preserves_records_and_lint(tmp_path):
    import gzip

    from vdjdb.qc.lint import lint_file

    plain = _chunk(tmp_path, 'PMID_1.tsv', [{'cdr3.beta': 'CASSLGQETQYF',
                                          'reference.id': 'PMID:1'}])
    compressed = tmp_path / 'PMID_1.tsv.gz'
    compressed.write_bytes(gzip.compress(plain.read_bytes(), mtime=0))
    assert read_chunk(compressed).drop('chunk.file').equals(read_chunk(plain).drop('chunk.file'))
    assert [(x.code, x.detail) for x in lint_file(compressed)] == [
        (x.code, x.detail) for x in lint_file(plain)]
    assert chunk_files(tmp_path) == [plain, compressed]
