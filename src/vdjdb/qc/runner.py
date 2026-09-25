"""Drive the chunk checks and report.

Phase 0 runs the text-level lint only. Phase 3 (``feature/io-qc``) adds the vectorised table rules
ported from ``py_src/ChunkQC.py``.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

from ..config import Paths
from .lint import Finding, lint


def _default_chunks() -> list[Path]:
    chunks = Paths.discover().chunks
    return sorted(p for p in chunks.iterdir() if p.suffix in {".txt", ".tsv"} and not p.name.startswith("."))


def run_qc(paths: list[Path] | None, *, strict: bool = True, report: Path | None = None) -> int:
    """Return an exit code: 0 clean, 1 findings under ``strict``."""
    targets = list(paths) if paths else _default_chunks()
    findings: list[Finding] = lint(targets)

    if report is not None:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(
            "file\tcode\tdetail\n" + "".join(f"{f}\n" for f in findings), encoding="utf-8"
        )

    by_code = Counter(f.code for f in findings)
    print(f"chunks checked: {len(targets)}")
    if not findings:
        print("no findings")
        return 0

    print(f"findings: {len(findings)}")
    for code, n in by_code.most_common():
        print(f"  {n:6d}  {code}")
    for f in findings[:40]:
        print(f"    {f}", file=sys.stderr)
    if len(findings) > 40:
        print(f"    ... and {len(findings) - 40} more", file=sys.stderr)

    return 1 if strict else 0
