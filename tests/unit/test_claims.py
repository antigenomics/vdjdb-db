"""Claims the repository makes about itself, checked instead of trusted.

Two defects found by accident on 2026-09-28 were the same shape: **a comment asserting a property the
code had stopped having.**

* `pyproject.toml` declared a pytest marker for "runtime and peak-RSS budgets; needs RUN_BENCHMARK=1".
  No test ever used it, and `RUN_BENCHMARK` appeared in exactly one place in the repository - that
  marker's own description.
* The same file said `validate/metrics_lib.py` is "vendored verbatim and must stay so", which stopped
  being true the moment #606 put a deliberate tie-break divergence in it.

Neither was caught by anything, because a comment cannot fail. Each test here turns one such claim
into something that can, chosen because the claim is load-bearing and the check is cheap - not to
police prose in general.
"""
from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "vdjdb"
TESTS = ROOT / "tests"


def _pyproject() -> dict:
    return tomllib.loads((ROOT / "pyproject.toml").read_bytes().decode())


def test_every_declared_pytest_marker_is_used_by_some_test() -> None:
    """A marker nobody applies is a promise nobody keeps -- `benchmark` promised a peak-RSS budget and
    had no test for the whole of phase 4."""
    declared = {m.split(":")[0].strip()
                for m in _pyproject()["tool"]["pytest"]["ini_options"]["markers"]}
    used = set()
    for p in TESTS.rglob("test_*.py"):
        used |= set(re.findall(r"mark\.(\w+)", p.read_text()))
    unused = declared - used
    assert not unused, (
        f"declared in pyproject.toml and applied by no test: {sorted(unused)}. Either write the test "
        f"or delete the marker; a marker that names a budget nobody measures reads as coverage.")


def test_pandas_stays_out_of_the_build() -> None:
    """Every stage from `chunks/` to the three shipped formats is polars.

    pandas is confined to `validate/`, where `metrics_lib.py` is vendored from TCREMP for benchmark
    comparability and one `to_pandas()` feeds it. It is in the `motifs` extra, not the core
    dependencies, and `pyproject.toml` says the build must not acquire one. This is that sentence, as
    a check.
    """
    offenders = sorted(
        str(p.relative_to(SRC)) for p in SRC.rglob("*.py")
        if p.parts[-2] != "validate"
        and re.search(r"^\s*import pandas|^\s*from pandas|\.to_pandas\(", p.read_text(), re.M))
    assert not offenders, (
        f"pandas outside validate/: {offenders}. The build is polars; see the note beside the "
        f"`motifs` extra in pyproject.toml.")


def test_pandas_is_not_a_core_dependency() -> None:
    core = " ".join(_pyproject()["project"]["dependencies"])
    assert "pandas" not in core, "pandas must stay in the `motifs` extra, never in the core deps"


def test_every_vdjdb_command_the_docs_name_exists() -> None:
    """Docs drift the moment a command is renamed, and a published command that does not run is worse
    than an undocumented one."""
    from vdjdb.cli import app

    known = {c.name or c.callback.__name__ for c in app.registered_commands}
    known |= {f"{g.name} {c.name}" for g in app.registered_groups
              for c in g.typer_instance.registered_commands}
    known |= {g.name for g in app.registered_groups}

    named: dict[str, set[str]] = {}
    for page in sorted((ROOT / "docs").rglob("*.md")):
        for m in re.findall(r"vdjdb ([a-z][a-z-]*(?: [a-z][a-z-]*)?)", page.read_text()):
            # Two words may be `group command` or `command <arg>`; accept either reading.
            if m in known or m.split()[0] in known:
                continue
            named.setdefault(m, set()).add(str(page.relative_to(ROOT)))
    unknown = {k: sorted(v) for k, v in named.items()
               if not any(k.startswith(c) for c in known)}
    assert not unknown, f"docs name commands the CLI does not have: {unknown}"


@pytest.mark.parametrize("rules_file", ["expected_diffs.toml", "motif_metrics.tsv",
                                        "build_timings.tsv", "motif_timings.tsv",
                                        "qc_advisories.tsv"])
def test_every_declared_baseline_exists_and_is_not_empty(rules_file: str) -> None:
    """`rules/` holds the declarations every gate compares against. A gate whose baseline vanished
    passes on an empty join, which is the failure mode that reads as success."""
    p = ROOT / "rules" / rules_file
    assert p.exists(), f"{p} is named by a gate and missing"
    assert p.stat().st_size > 0, f"{p} is empty; a gate with an empty baseline gates nothing"
