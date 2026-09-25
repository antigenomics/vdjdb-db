"""The difference ledger -- the instrument every later phase is measured with.

``vdjdb diff <reference> <candidate>`` compares a released bundle against a freshly built one in
three passes:

1. **file set** -- exact names and count;
2. **two digests per file** -- raw sha256 and *canonical* sha256 (data rows sorted). Canonical
   equality is the gate; raw equality is informational, because ``runBuidDatabase.py`` iterates
   ``os.listdir`` and readdir order is host-dependent, so the 2026-06-03 release cannot be
   reproduced byte-for-byte by anyone, including the pipeline that produced it;
3. **row-level classification** -- rows keyed on a stable identity, bucketed into
   ``only-in-reference`` / ``only-in-candidate`` / ``changed`` / ``identical``, and every changed
   **cell** attributed to a declared rule in ``rules/expected_diffs.toml``.

The gate is deliberately two-sided: an unattributed difference fails, and **a rule that fires a
different number of times than it declares also fails**. Without the second half a rule is a
description rather than a measurement, and a regression that happens to touch a declared column
would pass silently.
"""
from __future__ import annotations

import hashlib
import tomllib
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from ..schema import FULL_COLUMNS, SLIM_COLUMNS, VDJDB_COLUMNS

#: Row identity per tabular file. Keys are **not** unique -- one paper reporting the same TCR in
#: several donors is several rows -- so rows are compared as a multiset within each key group.
KEYS: dict[str, tuple[str, ...]] = {
    "vdjdb.txt": ("gene", "cdr3", "v.segm", "j.segm", "species",
                  "mhc.a", "mhc.b", "antigen.epitope", "reference.id"),
    "vdjdb.slim.txt": ("gene", "cdr3", "v.segm", "j.segm", "species",
                       "mhc.a", "mhc.b", "antigen.epitope", "reference.id"),
    "vdjdb_full.txt": ("cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta",
                       "species", "mhc.a", "mhc.b", "antigen.epitope", "reference.id"),
    "cluster_members.txt": ("cid", "cdr3aa", "gene", "species", "antigen.epitope"),
    "motif_pwms.txt": ("cid", "pos", "aa", "len", "gene", "species", "antigen.epitope"),
    "cluster_members_tcremp.txt": ("cid", "cdr3aa", "gene", "species", "antigen.epitope"),
    "motif_pwms_tcremp.txt": ("cid", "pos", "aa", "len", "gene", "species", "antigen.epitope"),
}

#: Expected column count, asserted before anything else. A width change means a positional-parse
#: break in ``vdjdb-web`` (CLAUDE.md, hard rule 1), not a row-level difference.
WIDTHS: dict[str, int] = {
    "vdjdb.txt": len(VDJDB_COLUMNS),
    "vdjdb.slim.txt": len(SLIM_COLUMNS),
    "vdjdb_full.txt": len(FULL_COLUMNS),
    "cluster_members.txt": 19,
    "motif_pwms.txt": 27,
}

#: Compared line-by-line rather than row-by-row: they are short, and their rows are the schema.
LINEWISE = ("vdjdb.meta.txt", "vdjdb.slim.meta.txt", "latest-version.txt")


# ---------------------------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class CellDiff:
    """One changed cell. ``key`` is the row identity, not a line number."""

    file: str
    column: str
    old: str
    new: str
    key: str


@dataclass(frozen=True, slots=True)
class FileReport:
    name: str
    raw_equal: bool
    canonical_equal: bool
    rows_reference: int = 0
    rows_candidate: int = 0
    only_in_reference: int = 0
    only_in_candidate: int = 0
    changed_rows: int = 0
    cells: tuple[CellDiff, ...] = ()
    note: str = ""


@dataclass
class DiffReport:
    missing: list[str] = field(default_factory=list)
    added: list[str] = field(default_factory=list)
    files: list[FileReport] = field(default_factory=list)
    unattributed: list[CellDiff] = field(default_factory=list)
    rule_counts: dict[str, int] = field(default_factory=dict)
    rule_expected: dict[str, int] = field(default_factory=dict)

    @property
    def miscounted_rules(self) -> list[str]:
        """Rules that fired a different number of times than they declare."""
        return sorted(r for r, want in self.rule_expected.items()
                      if want >= 0 and self.rule_counts.get(r, 0) != want)

    @property
    def ok(self) -> bool:
        return not (self.missing or self.added or self.unattributed or self.miscounted_rules
                    or any(not f.canonical_equal and not f.cells and f.name not in LINEWISE
                           for f in self.files))


# ---------------------------------------------------------------------------------------------
# Bundles
# ---------------------------------------------------------------------------------------------

class Bundle:
    """A release bundle: a zip (any single top-level directory is stripped) or a directory."""

    def __init__(self, path: Path) -> None:
        self.path = path
        if path.is_dir():
            self._zip = None
            self._members = {p.name: p for p in sorted(path.rglob("*")) if p.is_file()}
        else:
            self._zip = zipfile.ZipFile(path)
            self._members = {Path(n).name: n for n in self._zip.namelist()      # type: ignore[misc]
                             if not n.endswith("/")}

    @property
    def names(self) -> list[str]:
        return sorted(self._members)

    def read_bytes(self, name: str) -> bytes:
        m = self._members[name]
        return self._zip.read(m) if self._zip else Path(m).read_bytes()   # type: ignore[arg-type]


def _digests(data: bytes) -> tuple[str, str]:
    """``(raw, canonical)`` sha256. Canonical sorts the data rows, keeping the header first."""
    raw = hashlib.sha256(data).hexdigest()
    lines = data.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    head, body = (lines[:1], lines[1:]) if lines else ([], [])
    canonical = hashlib.sha256(b"\n".join(head + sorted(body))).hexdigest()
    return raw, canonical


def _read_table(data: bytes) -> pl.DataFrame:
    """Every cell as a string. No type inference: ``0`` and ``0.0`` are different cells here."""
    return pl.read_csv(data, separator="\t", quote_char=None, has_header=True,
                       infer_schema_length=0, truncate_ragged_lines=False).fill_null("")


# ---------------------------------------------------------------------------------------------
# Row-level comparison
# ---------------------------------------------------------------------------------------------

def _rows_by_key(df: pl.DataFrame, key: tuple[str, ...]) -> dict[str, list[tuple[str, ...]]]:
    present = [c for c in key if c in df.columns]
    keys = df.select(pl.concat_str(present, separator="\x1f")).to_series().to_list()
    out: dict[str, list[tuple[str, ...]]] = defaultdict(list)
    for k, row in zip(keys, df.iter_rows(), strict=True):
        out[k].append(row)
    return out


def _compare_table(name: str, ref: bytes, cand: bytes) -> FileReport:
    raw_r, can_r = _digests(ref)
    raw_c, can_c = _digests(cand)
    a, b = _read_table(ref), _read_table(cand)

    want = WIDTHS.get(name)
    if want is not None and len(a.columns) != want:
        # The reference itself is not the shape the registry declares -- fail loudly rather than
        # comparing against a file we have misidentified.
        return FileReport(name, raw_r == raw_c, False, a.height, b.height,
                          note=f"reference has {len(a.columns)} columns, registry declares {want}")

    if a.columns != b.columns:
        gone, new = set(a.columns) - set(b.columns), set(b.columns) - set(a.columns)
        note = (f"column set differs: -{sorted(gone)} +{sorted(new)}" if (gone or new)
                else f"columns reordered: {a.columns} -> {b.columns}")
        return FileReport(name, raw_r == raw_c, False, a.height, b.height, note=note)

    key = KEYS.get(name)
    if key is None:
        return FileReport(name, raw_r == raw_c, can_r == can_c, a.height, b.height,
                          note="no identity key declared; digest comparison only")

    ref_rows, cand_rows = _rows_by_key(a, key), _rows_by_key(b, key)
    cols = a.columns
    only_ref = only_cand = changed = 0
    cells: list[CellDiff] = []

    for k in set(ref_rows) | set(cand_rows):
        left, right = Counter(ref_rows.get(k, [])), Counter(cand_rows.get(k, []))
        common = left & right                      # identical rows, in multiset arithmetic
        left, right = sorted((left - common).elements()), sorted((right - common).elements())
        # Pair what is left positionally: same key, same count, differing content is a change.
        for lrow, rrow in zip(left, right, strict=False):
            changed += 1
            cells.extend(CellDiff(name, c, lv, rv, k)
                         for c, lv, rv in zip(cols, lrow, rrow, strict=True) if lv != rv)
        only_ref += max(len(left) - len(right), 0)
        only_cand += max(len(right) - len(left), 0)

    return FileReport(name, raw_r == raw_c, can_r == can_c, a.height, b.height,
                      only_ref, only_cand, changed, tuple(cells))


def _compare_lines(name: str, ref: bytes, cand: bytes) -> FileReport:
    raw_r, can_r = _digests(ref)
    raw_c, can_c = _digests(cand)
    a = ref.decode().splitlines()
    b = cand.decode().splitlines()
    cells = tuple(CellDiff(name, "line", old, new, str(i))
                  for i, (old, new) in enumerate(zip(a, b, strict=False)) if old != new)
    return FileReport(name, raw_r == raw_c, can_r == can_c, len(a), len(b),
                      max(len(a) - len(b), 0), max(len(b) - len(a), 0), len(cells), cells)


def _compare_opaque(name: str, ref: bytes, cand: bytes) -> FileReport:
    raw_r, can_r = _digests(ref)
    raw_c, can_c = _digests(cand)
    return FileReport(name, raw_r == raw_c, can_r == can_c,
                      note="" if can_r == can_c else "content differs (not a table)")


# ---------------------------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class Rule:
    id: str
    file: str | None = None
    column: str | None = None
    from_: str | None = None
    to: str | None = None
    rows: int = -1      # -1 = "no declared count"; see DiffReport.miscounted_rules

    def matches(self, c: CellDiff) -> bool:
        return ((self.file is None or self.file == c.file)
                and (self.column is None or self.column == c.column)
                and (self.from_ is None or self.from_ == c.old)
                and (self.to is None or self.to == c.new))


def load_rules(path: Path) -> list[Rule]:
    if not path.exists():
        return []
    raw = tomllib.loads(path.read_text())
    return [Rule(id=r["id"], file=r.get("file"), column=r.get("column"),
                 from_=r.get("from"), to=r.get("to"), rows=int(r.get("rows", -1)))
            for r in raw.get("rule", [])]


# ---------------------------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------------------------

def diff(reference: Path, candidate: Path, rules_path: Path | None = None) -> DiffReport:
    """Compare two bundles and attribute every changed cell to a declared rule."""
    ref, cand = Bundle(reference), Bundle(candidate)
    rules = load_rules(rules_path) if rules_path else []
    report = DiffReport(
        missing=[n for n in ref.names if n not in cand.names],
        added=[n for n in cand.names if n not in ref.names],
        rule_expected={r.id: r.rows for r in rules},
    )

    for name in ref.names:
        if name in report.missing:
            continue
        a, b = ref.read_bytes(name), cand.read_bytes(name)
        if name in KEYS or name in WIDTHS:
            report.files.append(_compare_table(name, a, b))
        elif name in LINEWISE:
            report.files.append(_compare_lines(name, a, b))
        else:
            report.files.append(_compare_opaque(name, a, b))

    for f in report.files:
        for cell in f.cells:
            for rule in rules:
                if rule.matches(cell):
                    report.rule_counts[rule.id] = report.rule_counts.get(rule.id, 0) + 1
                    break
            else:
                report.unattributed.append(cell)
    return report


def render(report: DiffReport) -> str:
    """The release-notes table. ``--report`` writes this, so the notes write themselves."""
    out = ["# Difference ledger", ""]
    if report.missing:
        out += [f"**Missing from the candidate:** {', '.join(report.missing)}", ""]
    if report.added:
        out += [f"**New in the candidate:** {', '.join(report.added)}", ""]

    out += ["| File | Raw | Canonical | Rows ref | Rows cand | Only ref | Only cand | Changed |",
            "|---|---|---|---|---|---|---|---|"]
    for f in sorted(report.files, key=lambda f: f.name):
        tick = lambda v: "yes" if v else "no"   # noqa: E731
        out.append(f"| `{f.name}` | {tick(f.raw_equal)} | {tick(f.canonical_equal)} | "
                   f"{f.rows_reference} | {f.rows_candidate} | {f.only_in_reference} | "
                   f"{f.only_in_candidate} | {f.changed_rows} |")
        if f.note:
            out.append(f"| | | | | | | | _{f.note}_ |")

    if report.rule_expected:
        out += ["", "## Declared rules", "", "| Rule | Declared | Fired |", "|---|---|---|"]
        for rid, want in sorted(report.rule_expected.items()):
            got = report.rule_counts.get(rid, 0)
            flag = "" if want < 0 or want == got else "  **MISCOUNT**"
            out.append(f"| `{rid}` | {'—' if want < 0 else want} | {got}{flag} |")

    if report.unattributed:
        by_col = Counter((c.file, c.column) for c in report.unattributed)
        out += ["", f"## Unattributed differences — {len(report.unattributed)} cells", "",
                "| File | Column | Cells | Example |", "|---|---|---|---|"]
        for (fname, col), n in by_col.most_common(40):
            ex = next(c for c in report.unattributed if c.file == fname and c.column == col)
            out.append(f"| `{fname}` | `{col}` | {n} | `{ex.old[:40]}` → `{ex.new[:40]}` |")

    out += ["", f"**Verdict: {'PASS' if report.ok else 'FAIL'}**"]
    return "\n".join(out) + "\n"
