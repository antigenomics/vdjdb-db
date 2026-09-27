"""Does this build still reproduce the last release?

Marked ``release``: needs ``VDJDB_REFERENCE_ZIP`` and a built legacy directory
(``VDJDB_LEGACY_DIR``, default ``out/legacy``).

``tests/release/test_reference_contract.py`` checks the *reference* -- its shape, its column orders,
its zero point against itself. This module runs the comparison that the reference exists for, and
asserts the three properties that make it a gate rather than a report:

1. every changed cell is attributed to a declared rule,
2. every rule fires exactly the number of rows it declares,
3. every row a file gained or lost is declared too.

The second is what turns a rule into a measurement. Without it a regression that happened to touch a
column some rule already covers would be absorbed in silence, and the run would still pass.

Why it is a test and not only the ``vdjdb diff`` step in ``build.yml``: the workflow runs the
comparison on the whole bundle and prints a verdict, which says *that* something is undeclared. The
assertions below say which property failed, and they run on a developer's machine against a build
they already have.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from vdjdb.compare.diff import Bundle, DiffReport, diff
from vdjdb.schema import SLIM_COLUMNS, VDJDB_COLUMNS

pytestmark = pytest.mark.release

#: The members the assembly stage produces. The motif files and the dashboard fragment come from
#: later stages, so a legacy build is a deliberately partial candidate.
LEGACY_MEMBERS = ("vdjdb.txt", "vdjdb.slim.txt", "vdjdb_full.txt",
                  "vdjdb.meta.txt", "vdjdb.slim.meta.txt")

DATA_TABLES = ("vdjdb.txt", "vdjdb.slim.txt", "vdjdb_full.txt")

#: name -> the table whose header the file's metadata must reproduce (hard rule 2).
METADATA_OF = {"vdjdb.meta.txt": ("vdjdb.txt", VDJDB_COLUMNS),
               "vdjdb.slim.meta.txt": ("vdjdb.slim.txt", SLIM_COLUMNS)}

RULES = Path("rules/expected_diffs.toml")


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
def candidate() -> Path:
    d = Path(os.environ.get("VDJDB_LEGACY_DIR", "out/legacy"))
    if not (d / "vdjdb.txt").exists():
        pytest.skip(f"no legacy build at {d}; run `vdjdb build --out out/`")
    return d


@pytest.fixture(scope="module")
def report(reference: Path, candidate: Path) -> DiffReport:
    return diff(reference, candidate, RULES, only=LEGACY_MEMBERS)


# --------------------------------------------------------------------------------------------
# The gate
# --------------------------------------------------------------------------------------------

def test_every_changed_cell_is_attributed(report: DiffReport) -> None:
    detail = "\n".join(f"  {c.file} {c.column} [{c.key}]: {c.old!r} -> {c.new!r}"
                       for c in report.unattributed[:20])
    assert not report.unattributed, \
        f"{len(report.unattributed)} undeclared cell differences:\n{detail}"


def test_every_rule_fires_the_number_of_rows_it_declares(report: DiffReport) -> None:
    detail = "\n".join(f"  {r}: declared {report.rule_expected[r]}, fired "
                       f"{report.rule_counts.get(r, 0)}" for r in report.miscounted_rules)
    assert not report.miscounted_rules, f"miscounted rules:\n{detail}"


def test_no_declaration_outlives_the_data_it_describes(report: DiffReport) -> None:
    """A rename that matches nothing is a rule about a release that is no longer the reference."""
    assert not report.stale_renames


def test_every_row_gained_or_lost_is_declared(report: DiffReport) -> None:
    for f in report.files:
        d = report.row_deltas.get(f.name)
        want = (d.removed, d.added) if d else (0, 0)
        assert (f.only_in_reference, f.only_in_candidate) == want, (
            f"{f.name}: {f.only_in_reference} rows only in the reference and "
            f"{f.only_in_candidate} only in the candidate; declared {want}")


def test_the_comparison_passes(report: DiffReport) -> None:
    """The one assertion the release job makes. The four above say which way it failed."""
    assert report.ok


# --------------------------------------------------------------------------------------------
# That the comparison measured what it claims to
# --------------------------------------------------------------------------------------------

def test_the_bundle_contains_every_legacy_member(report: DiffReport) -> None:
    assert not report.missing, f"the candidate is missing {report.missing}"
    assert not report.added, f"the candidate ships {report.added}, which the release does not"


def test_every_member_was_compared_row_by_row(report: DiffReport) -> None:
    """A file compared on its digest alone is a file whose differences were never enumerated.

    `_compare_table` falls back to a digest when a width or a column set disagrees, and reports it
    in `note`. A pass built on five digest comparisons would mean nothing.
    """
    for f in report.files:
        assert f.compared_rows, f"{f.name}: {f.note or 'no row comparison ran'}"
    assert {f.name for f in report.files} == set(LEGACY_MEMBERS)


@pytest.mark.parametrize("name", DATA_TABLES)
def test_the_header_is_unchanged_from_the_release(name: str, reference: Path,
                                                  candidate: Path) -> None:
    """`vdjdb-web` reads these positionally, and the release is what its parser was built against."""
    ref = Bundle(reference).read_bytes(name).split(b"\n", 1)[0].decode().split("\t")
    cand = Bundle(candidate).read_bytes(name).split(b"\n", 1)[0].decode().split("\t")
    assert cand == ref


@pytest.mark.parametrize("meta", sorted(METADATA_OF))
def test_the_metadata_describes_the_table_that_shipped_beside_it(meta: str,
                                                                 candidate: Path) -> None:
    """Hard rule 2, asserted on the build rather than on the registry.

    ``test_schema.py`` ties the metadata to the header the registry declares. This ties it to the
    bytes that would go in the zip, which is the file ``vdjdb-web`` builds its schema from.
    """
    table, declared = METADATA_OF[meta]
    b = Bundle(candidate)
    names = [line.split("\t")[0] for line in b.read_bytes(meta).decode().splitlines()[1:]]
    header = b.read_bytes(table).split(b"\n", 1)[0].decode().split("\t")
    assert names == header
    assert names == list(declared)
