"""Drive the chunk checks and report.

Two tiers, in order of cost:

1. :mod:`vdjdb.qc.lint` -- text-level, never parses the body. Catches the class of problem that
   otherwise appears as a confusing parser error three stages later.
2. :mod:`vdjdb.qc.rules` -- row-level, vectorised over all 230 chunks in one pass.

Lint findings are warnings: 99 of the 230 chunks are CRLF and two have malformed headers, and a hard
gate on those would block every unrelated submission until the `.tsv` migration (#497) lands. Rule
findings are errors and fail under ``--strict``, restoring the behaviour the Groovy build had and the
Python port replaced with ``warnings.warn``. Measured: 230 of 230 chunks pass today, so there is no
quarantine list to grandfather.

``duplicate`` is neither: it is a curation signal. Per-chunk deduplication removes those rows on
the way in, and their count is the gap between the 203,308 raw rows and the released 192,753.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import polars as pl

from ..config import Paths
from ..io.chunks import chunk_files, read_chunks
from .lint import Finding, lint
from .rules import check, summarise

#: Reported, never fatal. See the module docstring.
ADVISORY = frozenset({"crlf", "unknown-column", "reference-id-form", "duplicate",
                      # #561: only a curator can decide which of the two chains is the wrong one.
                      "alpha and beta cdr3 identical"})


def _report_frame(lint_findings: list[Finding], row_findings: pl.DataFrame) -> pl.DataFrame:
    text = pl.DataFrame(
        {"file": [f.file for f in lint_findings],
         "row": [0] * len(lint_findings),
         "code": [f.code for f in lint_findings],
         "detail": [f.detail for f in lint_findings]},
        schema={"file": pl.Utf8, "row": pl.UInt32, "code": pl.Utf8, "detail": pl.Utf8},
    )
    rows = row_findings.select(
        pl.col("chunk.file").alias("file"),
        pl.col("chunk.row").alias("row"),
        pl.col("rule").alias("code"),
        pl.lit("").alias("detail"),
    )
    return pl.concat([text, rows]).sort("file", "row", "code")


def run_qc(paths: list[Path] | None, *, strict: bool = True, report: Path | None = None) -> int:
    """Return an exit code: 0 clean, 1 if any fatal finding was raised under ``strict``."""
    targets = list(paths) if paths else chunk_files(Paths.discover().chunks)
    lint_findings = lint(targets)
    try:
        rows = read_chunks(targets, deduplicate=False)
        row_findings = check(rows)
    except Exception as exc:
        # A chunk whose header is wrong cannot be parsed as a table at all, and that is the most
        # common submission mistake. Raising here produced a traceback and skipped the `--report`
        # write below, so `chunk-check` uploaded no report and its pull-request comment said
        # "No QC findings" for the case that most needed a finding.
        lint_findings = [*lint_findings, Finding(file=targets[0].name if targets else "",
                                                 code="unreadable", detail=str(exc).strip())]
        rows = pl.DataFrame()
        row_findings = pl.DataFrame(schema={"chunk.file": pl.Utf8, "chunk.row": pl.UInt32,
                                            "rule": pl.Utf8})

    if report is not None:
        report.parent.mkdir(parents=True, exist_ok=True)
        _report_frame(lint_findings, row_findings).write_csv(report, separator="\t")

    print(f"chunks checked: {len(targets)}    rows: {rows.height}")

    by_code = Counter(f.code for f in lint_findings)
    if by_code:
        print("text-level findings:")
        for code, n in by_code.most_common():
            print(f"  {n:6d}  {code}{'' if code not in ADVISORY else '   (advisory)'}")

    if row_findings.height:
        print("row-level findings:")
        for row in summarise(row_findings).iter_rows(named=True):
            advisory = "   (advisory)" if row["rule"] in ADVISORY else ""
            print(f"  {row['rows']:6d}  {row['rule']}   in {row['chunks']} chunk(s){advisory}")
    else:
        print("row-level findings: none")

    fatal = (sum(n for c, n in by_code.items() if c not in ADVISORY)
             + row_findings.filter(~pl.col("rule").is_in(list(ADVISORY))).height)
    if fatal and strict:
        print(f"FAILED: {fatal} fatal finding(s)")
        return 1
    if fatal:
        print(f"{fatal} fatal finding(s), not failing because --no-strict")
    return 0
