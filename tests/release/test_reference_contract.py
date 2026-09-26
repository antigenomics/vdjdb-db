"""The registry against a real VDJdb release zip.

Marked ``release``: needs ``VDJDB_REFERENCE_ZIP`` and runs nightly, not per-PR. It is the check that
the positional contracts in ``vdjdb.schema`` describe the file consumers actually parse, rather than
what the source code says they should.
"""
from __future__ import annotations

import os
import zipfile
from pathlib import Path

import pytest

from vdjdb.compare.diff import Bundle, diff
from vdjdb.schema import TABLES

pytestmark = pytest.mark.release

FILE_TO_TABLE = {
    "vdjdb.txt": "vdjdb",
    "vdjdb.slim.txt": "slim",
    "vdjdb_full.txt": "full",
    "cluster_members.txt": "cluster_members",
    "motif_pwms.txt": "motif_pwms",
}

#: Row counts of the 2026-06-03 release. A change here is a real change in the database, not a bug.
EXPECTED_ROWS = {
    "vdjdb.txt": 284_546,
    "vdjdb.slim.txt": 197_729,
    "vdjdb_full.txt": 192_753,
    "cluster_members.txt": 55_636,
    "motif_pwms.txt": 40_061,
}


@pytest.fixture(scope="module")
def reference() -> Path:
    env = os.environ.get("VDJDB_REFERENCE_ZIP")
    if not env:
        pytest.skip("set VDJDB_REFERENCE_ZIP to a release zip or extracted directory")
    p = Path(env)
    if not p.exists():
        pytest.skip(f"{p} does not exist")
    return p


@pytest.fixture(scope="module")
def headers(reference: Path) -> dict[str, list[str]]:
    b = Bundle(reference)
    return {n: b.read_bytes(n).split(b"\n", 1)[0].decode().split("\t")
            for n in FILE_TO_TABLE if n in b.names}


@pytest.mark.parametrize("name", sorted(FILE_TO_TABLE))
def test_column_order_matches_the_registry(name: str, headers: dict[str, list[str]]) -> None:
    """``vdjdb-web`` parses the motif files positionally with no header check (ROADMAP.md rule 1)."""
    if name not in headers:
        pytest.skip(f"{name} not in this bundle")
    assert headers[name] == list(TABLES[FILE_TO_TABLE[name]])


def test_the_bundle_has_exactly_the_specified_members(reference: Path) -> None:
    """``docs/outputs.md`` section 2: ten members, and the basenames are a vdjmatch contract."""
    names = set(Bundle(reference).names)
    required = {"vdjdb.txt", "vdjdb.meta.txt", "vdjdb.slim.txt", "vdjdb.slim.meta.txt",
                "vdjdb_full.txt", "cluster_members.txt", "motif_pwms.txt",
                "vdjdb_summary_embed.html", "LICENSE", "latest-version.txt"}
    assert required <= names, f"missing: {sorted(required - names)}"


@pytest.mark.parametrize("name", ["vdjdb.txt", "vdjdb.slim.txt", "vdjdb_full.txt"])
def test_no_field_is_csv_quoted(name: str, reference: Path) -> None:
    """Why ``quote_style="never"`` reproduces the release exactly.

    The files *do* contain ``"`` characters -- every ``method`` / ``meta`` / ``cdr3fix`` cell is
    JSON. What they never contain is a *quoted field*: a default CSV writer would wrap those cells
    in quotes and double the inner ones, changing every JSON cell in the database. So the property
    is that no field begins with a quote, not that no quote exists.
    """
    for line in Bundle(reference).read_bytes(name).decode().split("\n"):
        if not line:
            continue
        assert not any(f.startswith('"') for f in line.split("\t")), line[:120]


def test_vdjdb_txt_is_rectangular_with_no_empty_cdr3fix(reference: Path) -> None:
    """``vdjdb-web`` right-pads short lines and crashes on ``Json.parse("")``; both must stay moot."""
    data = Bundle(reference).read_bytes("vdjdb.txt").decode().split("\n")
    header = data[0].split("\t")
    fix = header.index("cdr3fix")
    ragged = empty = 0
    for line in data[1:]:
        if not line:
            continue
        f = line.split("\t")
        ragged += len(f) != len(header)
        empty += not f[fix]
    assert ragged == 0
    assert empty == 0


@pytest.mark.parametrize("name", sorted(EXPECTED_ROWS))
def test_row_counts(name: str, reference: Path) -> None:
    b = Bundle(reference)
    if name not in b.names:
        pytest.skip(f"{name} not in this bundle")
    data = b.read_bytes(name)
    rows = data.count(b"\n") - (0 if data.endswith(b"\n") else -1) - 1
    assert rows == EXPECTED_ROWS[name]


def test_the_ledger_finds_no_difference_between_a_release_and_itself(reference: Path) -> None:
    """The instrument's own zero point. If this drifts, nothing measured with it means anything."""
    report = diff(reference, reference)
    assert report.ok
    assert not report.unattributed
    assert all(f.raw_equal and f.canonical_equal for f in report.files)


def test_zip_members_live_under_a_single_directory(reference: Path) -> None:
    if reference.is_dir():
        pytest.skip("directory bundle")
    with zipfile.ZipFile(reference) as z:
        prefixes = {n.split("/")[0] for n in z.namelist()}
    assert len(prefixes) == 1, f"bundle has several top-level entries: {sorted(prefixes)}"
