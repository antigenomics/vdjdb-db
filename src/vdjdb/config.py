"""Paths, resolved once from the repo root.

Every legacy script in ``py_src/`` hard-codes ``../chunks`` and must be run from inside ``py_src/``.
Nothing here uses a relative path; the root is found by walking up for a marker file, so the CLI
works from any working directory.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_MARKERS = ("pyproject.toml", ".git")


def repo_root(start: Path | None = None) -> Path:
    """Walk up from ``start`` (default: this file) until a repo marker is found."""
    here = (start or Path(__file__)).resolve()
    for candidate in (here, *here.parents):
        if any((candidate / m).exists() for m in _MARKERS):
            return candidate
    raise RuntimeError(f"no repo root above {here}")


@dataclass(frozen=True, slots=True)
class Paths:
    """Every input and output location, derived from one root."""

    root: Path

    # inputs
    @property
    def chunks(self) -> Path:
        return self.root / "chunks"

    @property
    def patches(self) -> Path:
        return self.root / "patches"

    @property
    def res(self) -> Path:
        return self.root / "res"

    @property
    def proofreading(self) -> Path:
        return self.root / "proofreading"

    @property
    def summary(self) -> Path:
        return self.root / "summary"

    # outputs -- `out/`, not `build/`: `build/` is already gitignored as a Python packaging
    # convention and reusing it for release artifacts is confusing. See CLAUDE.md.
    @property
    def out(self) -> Path:
        return self.root / "out"

    @property
    def reports(self) -> Path:
        return self.out / "reports"

    @classmethod
    def discover(cls) -> Paths:
        return cls(repo_root())


#: The one seed for the whole build. Every stage that samples, shuffles, initialises a clustering
#: or hashes with a seed takes it from here -- never ``random.seed()`` at module scope, never an
#: unseeded default, never a per-call literal.
#:
#: A build that is merely *repeatable* (same answer twice on one host) is not enough: the difference
#: ledger declares rules with measured row counts, and a count is meaningless against a measurement
#: that moves between hosts or processes. Measured cost of getting this wrong: an unstable sort in
#: the ledger reported 158, 152 and 158 changed rows across three runs of the same comparison, and
#: manufactured five times the real difference.
SEED: int = 20260925
