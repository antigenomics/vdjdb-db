"""`ruff` over the tree, so `uv run pytest -q` predicts CI.

CI's `check` job runs `ruff check .` and the test suite did not, so a lint error was discoverable
only from a pull request: the `harmonise_vocabulary` import went in unsorted and an `E731` lambda
with it, both green locally across 1,149 tests and both red on the runner twenty minutes later.
That is the failure mode `CLAUDE.md` names for every other gate - a bar you cannot run locally gets
found at the worst moment - and ruff takes milliseconds, so there is no reason for it to be the
exception.

Skipped rather than failed when ruff is absent, because it lives in the `dev` extra and a `test`-only
install is a legitimate way to run this suite. CI installs both, so the gate is never silently off
where it matters.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def _ruff() -> list[str] | None:
    if shutil.which("ruff"):
        return ["ruff"]
    # Installed as a library into this interpreter's environment but not on PATH.
    probe = subprocess.run([sys.executable, "-m", "ruff", "--version"],
                           capture_output=True, cwd=ROOT, check=False)
    return [sys.executable, "-m", "ruff"] if probe.returncode == 0 else None


def test_ruff_check_passes() -> None:
    cmd = _ruff()
    if cmd is None:
        pytest.skip("ruff is not installed; it is in the `dev` extra and CI installs it")
    got = subprocess.run([*cmd, "check", "."], capture_output=True, text=True,
                         cwd=ROOT, check=False)
    assert got.returncode == 0, f"ruff check failed:\n{got.stdout}\n{got.stderr}"
