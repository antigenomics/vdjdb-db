"""Chunk lint: the cheap text-level checks, before anything is parsed as a table.

These run in milliseconds over all 230 chunks and catch the class of problem that otherwise appears
as a confusing pandas error three stages later: a file that is not UTF-8, a CRLF file, a header with
an empty or duplicated column name, a prose sentence used as a column name.

Measured on the corpus at the time of writing: 19 distinct header rows across 230 files, 99 files
CRLF, and four files with malformed headers. The fixes land with the `.tsv` rename (#497); until
then these are reported, and `--strict` decides whether they fail the build.
"""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

from ..schema import ALL_COLUMNS, COMPLEX_COLUMNS, KEPT_CURATION_COLUMNS

_BOM = b"\xef\xbb\xbf"
#: A column name longer than this is prose, not a name.
_MAX_NAME_LEN = 40

#: Accepted `reference.id` forms. Anything else is flagged for #347.
_REF_PREFIXES = ("PMID:", "doi:", "https://doi.org/", "http://", "https://")


@dataclass(frozen=True, slots=True)
class Finding:
    file: str
    code: str
    detail: str

    def __str__(self) -> str:
        return f"{self.file}\t{self.code}\t{self.detail}"


def lint_file(path: Path) -> list[Finding]:
    """Text-level checks on one chunk. Never parses the body as a table."""
    out: list[Finding] = []
    name = path.name
    raw = path.read_bytes()
    if path.suffix == ".gz":
        raw = gzip.decompress(raw)

    if raw.startswith(_BOM):
        out.append(Finding(name, "bom", "file starts with a UTF-8 BOM"))
        raw = raw[len(_BOM):]

    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        out.append(Finding(name, "encoding", f"not valid UTF-8: {exc}"))
        text = raw.decode("utf-8", errors="replace")

    if b"\r\n" in raw:
        out.append(Finding(name, "crlf", "CRLF line endings"))

    lines = text.splitlines()
    if not lines:
        out.append(Finding(name, "empty", "file has no lines"))
        return out

    header = lines[0].rstrip("\r").split("\t")

    for i, col in enumerate(header):
        if col == "":
            out.append(Finding(name, "empty-column-name", f"column {i} has an empty name"))
        elif len(col) > _MAX_NAME_LEN:
            out.append(Finding(name, "prose-column-name", f"column {i}: {col[:60]!r}..."))

    # `ChunkQC.check_exist` raised `ValueError("Empty file")` on a frame with no rows, which a
    # header-only submission produces. Nothing here reported it: `empty` needs a file with no lines at
    # all, and the reader returns a 0-row frame without complaint, so a chunk that contributes nothing
    # passed every gate.
    if len(lines) == 1:
        out.append(Finding(name, "no-data-rows", "header only, no records"))

    seen: set[str] = set()
    for col in header:
        if col in seen:
            out.append(Finding(name, "duplicate-column-name", col))
        seen.add(col)

    known = set(ALL_COLUMNS) | set(KEPT_CURATION_COLUMNS)
    for col in header:
        if col and col not in known and len(col) <= _MAX_NAME_LEN:
            out.append(Finding(name, "unknown-column", col))

    missing = [c for c in COMPLEX_COLUMNS if c not in seen]
    if missing:
        out.append(Finding(name, "missing-required-column", ",".join(missing)))

    if "reference.id" in header:
        idx = header.index("reference.id")
        bad: set[str] = set()
        for line in lines[1:]:
            fields = line.rstrip("\r").split("\t")
            if len(fields) > idx:
                ref = fields[idx].strip()
                if ref and not ref.startswith(_REF_PREFIXES):
                    bad.add(ref)
        for ref in sorted(bad)[:5]:
            out.append(Finding(name, "reference-id-form", ref))

    return out


def lint(paths: list[Path]) -> list[Finding]:
    """Text-level checks over several chunks.

    A path that does not exist is reported, not raised on. A pull request withdrawing a chunk names
    a file that is gone by the time the check runs, and `read_bytes` turned that into a
    `FileNotFoundError` traceback with no finding and no usable message.
    """
    findings: list[Finding] = []
    for p in sorted(paths):
        if not p.is_file():
            findings.append(Finding(p.name, "missing-file", f"{p} does not exist"))
            continue
        findings.extend(lint_file(p))
    return findings
