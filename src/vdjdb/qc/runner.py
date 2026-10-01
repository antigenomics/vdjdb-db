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
#: Findings that are reported but do not fail the run.
#:
#: The test for membership is whether the reader can produce the right record anyway. `chunks/` is
#: the submitters' data: the reader adapts to the shape a file arrives in and says what it found,
#: and a file is never edited into a shape the reader finds convenient. Rewriting the corpus to
#: silence a lint costs the thing the chunk-change rule exists to protect -- once every line of a
#: file has changed, a curation edit and a line-ending change are indistinguishable in `git log`.
#:
#: `prose-column-name` and `empty-column-name` are here for that reason, verified rather than
#: assumed. `PMID_24512815.tsv` carries two sentences of documentation as column names and
#: `PMID_40694338.tsv` opens with an unnamed column holding a row serial; columns are selected by
#: name, so both are ignored, and both files read with every field in the right column
#: (`cdr3.beta`, `v.beta`, `species`, `antigen.epitope` and `reference.id` all check out). A header
#: the reader cannot map is a different matter and still fails.
ADVISORY = frozenset({"crlf", "unknown-column", "reference-id-form", "duplicate",
                      # A named path that is gone: the caller passed a stale list, which is worth
                      # reporting but is not a defect in anybody's data.
                      "missing-file",
                      # The reader ignores columns it cannot name and reads the rest correctly.
                      "prose-column-name", "empty-column-name",
                      # #634. IMGT's F / ORF / P verdict on a named segment. A pseudogene call is not
                      # automatically wrong - a P gene can rearrange - and IMGT reclassifies genes
                      # between releases, so a gate here would fail on a reference update.
                      *(f"non-functional {c}" for c in
                        ("v.alpha", "j.alpha", "v.beta", "j.beta")),
                      # #561: only a curator can decide which of the two chains is the wrong one.
                      "alpha and beta cdr3 identical",
                      # #694, #625: a spreadsheet counter in a column that describes the antigen.
                      # Advisory because the repair is a patch entry or a chunk edit and either is a
                      # curation decision - and because the one current finding is already declared
                      # in `patches/mhc.dict`, so the rule is a regression guard rather than a gate
                      # on live data.
                      *(f"counter in {c}" for c in
                        ("antigen.gene", "antigen.species", "mhc.a", "mhc.b")),
                      # A named V or J whose chain has no CDR3. The call is information and the row
                      # is kept; the chain cannot reach an output, which is what this tells the
                      # submitter while they can still supply the sequence.
                      "segment call with no cdr3",
                      # A Cys after the one the junction opens with. Kept and flagged rather than
                      # gated: the Jurkat receptor carries one, so it is rare and not impossible, and
                      # 4,325 corpus rows over 108 chunks would fail a gate for a defect a curator has
                      # to read the source to settle.
                      "internal cysteine in cdr3.alpha", "internal cysteine in cdr3.beta",
                      # #402: a `meta.structure.id` that is not a PDB entry id. 2,765 rows in the
                      # corpus hold a figure or table reference there, and the field awards the top
                      # confidence score, so this has to be visible - but blanking them moves 6,004
                      # scores, which is a curation decision and not a submission error.
                      "structure id is not a PDB id",
                      # #637: a `method.identification` token `proofreading/method_vocabulary.tsv`
                      # has not settled. Advisory, and it must stay advisory - a submission naming a
                      # method nobody has seen is a method nobody has seen, not a defect. The corpus
                      # already carries `T-Scan`, `YAMTAD system` and `phage display`, none of which
                      # the specification page ever named. What the finding buys is that the next one
                      # is seen when it arrives.
                      "undeclared method.identification token",
                      # #696: a submitted `method.frequency` that its own count and total
                      # contradict. Which of the three the paper supports is a curation question and
                      # the repair is a chunk edit, so this reports. Zero findings today - the
                      # columns are new, so no submission has used them yet - which makes it a gate
                      # on the next one rather than a report on the corpus.
                      "frequency disagrees with its count and total"})


class _NothingToRead(Exception):
    """Every named path was missing. `lint` said so already; there is no table to check."""


def _summary_frame(lint_findings, row_findings: pl.DataFrame) -> pl.DataFrame:
    """One row per rule that fired: its level, how many findings and how many chunks.

    Sorted by ``(level, rule)`` rather than by count, so a diff between two builds lines up
    rule-for-rule instead of reordering when a count changes (hard rule 7 -- sort after anything
    unordered, and ``Counter.most_common`` is ordered by a value that moves).
    """
    rows = [{"level": "text", "rule": code, "findings": n, "chunks": -1,
             "advisory": code in ADVISORY}
            for code, n in Counter(f.code for f in lint_findings).items()]
    if row_findings.height:
        rows += [{"level": "row", "rule": r["rule"], "findings": r["rows"], "chunks": r["chunks"],
                  "advisory": r["rule"] in ADVISORY}
                 for r in summarise(row_findings).iter_rows(named=True)]
    schema = {"level": pl.Utf8, "rule": pl.Utf8, "findings": pl.Int64, "chunks": pl.Int64,
              "advisory": pl.Boolean}
    return pl.DataFrame(rows, schema=schema).sort("level", "rule")


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
    # `lint` has already reported any path that is gone. Drop them before reading, or a stale path
    # in the caller's list is diagnosed as `unreadable`, which means "this header cannot be parsed"
    # and sends a curator looking at a file that is not the problem.
    targets = [p for p in targets if p.is_file()]
    try:
        if not targets:
            raise _NothingToRead
        rows = read_chunks(targets, deduplicate=False)
        row_findings = check(rows)
    except _NothingToRead:
        rows = pl.DataFrame()
        row_findings = pl.DataFrame(schema={"chunk.file": pl.Utf8, "chunk.row": pl.UInt32,
                                            "rule": pl.Utf8})
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
        # The per-rule counts beside the per-row report. They were printed and nowhere else: an
        # advisory that jumped from 209 rows to 2,000 on a chunk merge would pass every gate and show
        # up only in a log line nobody diffs. `chunks/` is the data, so this is the curation signal.
        _summary_frame(lint_findings, row_findings).write_csv(
            report.with_name("qc-summary.tsv"), separator="\t")

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
