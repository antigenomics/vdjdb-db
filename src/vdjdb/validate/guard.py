"""Keep proprietary held-out data out of the repository and out of every release.

TCRvdb / MATCHMAKERS (Messemaker et al., doi:10.1101/2025.04.28.651095) is licensed for academic,
non-commercial use with **no redistribution, in whole or in part**. VDJdb ships AGPL-3.0-only. So
TCRvdb may be read during validation and must never be written anywhere this repository publishes.

It is also *held out*: motif clustering is tuned on the independent-study support count, and TCRvdb
is touched once at the end to validate. Tuning on it would invalidate it. See ROADMAP section 8 and
docs/outputs.md section 7.

This module is the enforcement, not the reminder:

- :func:`tcrvdb_path` is the only supported way to reach the file, and it reads an environment
  variable -- there is no default path inside the repository, so the file cannot be picked up by
  accident.
- :func:`scan_paths` fails a build if anything matching a proprietary fingerprint appears in the
  repo or in a bundle. CI runs it over the working tree and over every zip before publishing.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

ENV_VAR = "VDJDB_TCRVDB"

#: Filename fingerprints of the proprietary distribution. Matched case-insensitively on the whole
#: path, so a copy renamed into a subdirectory is still caught.
PROPRIETARY_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"tcrvdb", re.I),
    re.compile(r"matchmaker", re.I),
    re.compile(r"\d{2}_\d{2}_\d{4}_TCRvdb", re.I),
)

#: Column names that only appear together in the TCRvdb distribution. A file carrying all of these
#: is a copy of it whatever it is called.
_FINGERPRINT_COLUMNS = frozenset({"clonotype_aa", "epitope_aa", "hla_short", "padj", "baseMean"})

#: Files we read to check, but never large binaries.
_TEXT_SUFFIXES = frozenset({".tsv", ".csv", ".txt", ".json", ".md", ".parquet"})


@dataclass(frozen=True, slots=True)
class ProprietaryLeak:
    path: str
    reason: str

    def __str__(self) -> str:
        return f"{self.path}: {self.reason}"


def tcrvdb_path() -> Path:
    """Where the held-out TCRvdb table lives, from ``$VDJDB_TCRVDB``.

    Raises if unset or missing, rather than falling back to a path inside the repository -- a
    default here is exactly how a proprietary file ends up committed.
    """
    raw = os.environ.get(ENV_VAR)
    if not raw:
        raise RuntimeError(
            f"{ENV_VAR} is not set. TCRvdb is proprietary held-out validation data and is never "
            f"stored in this repository; point {ENV_VAR} at your own copy."
        )
    p = Path(raw).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"{ENV_VAR}={raw} does not exist")
    return p


def _header_is_tcrvdb(path: Path) -> bool:
    if path.suffix.lower() not in {".tsv", ".csv", ".txt"}:
        return False
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            header = fh.readline()
    except OSError:
        return False
    fields = set(re.split(r"[,\t]", header.strip()))
    return fields >= _FINGERPRINT_COLUMNS


def scan_paths(paths: list[Path], *, root: Path | None = None) -> list[ProprietaryLeak]:
    """Return every proprietary leak among ``paths``. Empty means clean."""
    leaks: list[ProprietaryLeak] = []
    for p in paths:
        rel = str(p.relative_to(root)) if root and p.is_relative_to(root) else str(p)
        for pat in PROPRIETARY_PATTERNS:
            if pat.search(rel):
                leaks.append(ProprietaryLeak(rel, f"filename matches proprietary pattern /{pat.pattern}/"))
                break
        else:
            if p.is_file() and p.suffix.lower() in _TEXT_SUFFIXES and _header_is_tcrvdb(p):
                leaks.append(ProprietaryLeak(rel, "header carries the TCRvdb column fingerprint"))
    return leaks
