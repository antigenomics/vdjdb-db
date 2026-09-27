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
    """The input directories, derived from one root.

    Only the three below: an output path is a CLI argument, so a property for it would be a
    fourth spelling of `--out`.
    """

    root: Path

    @property
    def chunks(self) -> Path:
        return self.root / "chunks"

    @property
    def patches(self) -> Path:
        return self.root / "patches"

    @property
    def res(self) -> Path:
        return self.root / "res"

    @classmethod
    def discover(cls) -> Paths:
        return cls(repo_root())


#: The one seed for the build. Every stage that samples, shuffles, initialises a clustering
#: or hashes with a seed takes it from here -- never ``random.seed()`` at module scope, never an
#: unseeded default, never a per-call literal.
#:
#: A build that is merely repeatable (same answer twice on one host) is not enough:
#: ``rules/expected_diffs.toml`` declares rules with measured row counts, and a count is meaningless
#: against a measurement that moves between hosts or processes. Measured cost of getting this wrong:
#: an unstable sort in ``vdjdb diff`` reported 158, 152 and 158 changed rows across three runs of the
#: same comparison, five times the number of differences actually present.
SEED: int = 20260925
