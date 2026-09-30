"""Comparison of a candidate build against a released bundle.

``vdjdb diff <reference> <candidate>`` compares a released bundle against a freshly built one in
three passes:

1. **file set** -- exact names and count;
2. **two digests per file** -- raw sha256 and *canonical* sha256 (data rows sorted). Canonical
   equality is the gate; raw equality is informational, because ``runBuidDatabase.py`` iterates
   ``os.listdir`` and readdir order is host-dependent, so the 2026-06-03 release cannot be
   reproduced byte-for-byte by anyone, including the pipeline that produced it;
3. **row-level classification** -- rows keyed on a stable identity, bucketed into
   ``only-in-reference`` / ``only-in-candidate`` / ``changed`` / ``identical``, and every changed
   cell attributed to a declared rule in ``rules/expected_diffs.toml``.

The gate is two-sided: an unattributed difference fails, and a rule that fires a different number
of times than it declares also fails. Without the second half a rule is a description rather than a
measurement, and a regression that touches a declared column would pass with no error.
"""
from __future__ import annotations

import hashlib
import json
import tomllib
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from ..schema import FULL_COLUMNS, SLIM_COLUMNS, VDJDB_COLUMNS

#: Row identity per tabular file. Keys are not unique -- one paper reporting the same TCR in
#: several donors is several rows -- so rows are compared as a multiset within each key group.
#:
#: The key must mirror the record identity, not merely the receptor. A narrower key was tried first
#: and produced 15,561 spurious cell differences against a build whose ``vdjdb_full.txt`` was
#: canonically identical: rows sharing a key were paired arbitrarily, so every field that
#: distinguishes two donors of one study reported as changed. Every one of those fields is in
#: :data:`CHUNK_DEDUP_KEY`, which is why the key now mirrors it.
KEYS: dict[str, tuple[str, ...]] = {
    # `antigen.gene` and `antigen.species` are deliberately absent: both are derived from
    # `antigen.epitope` by the patch table, so they annotate the epitope rather than identify the
    # record, and a nomenclature correction to either must read as a changed cell rather than as a
    # record removed and another added.
    "vdjdb.txt": ("gene", "cdr3", "v.segm", "j.segm", "species",
                  "mhc.a", "mhc.b", "mhc.class", "antigen.epitope", "reference.id"),
    "vdjdb.slim.txt": ("gene", "cdr3", "v.segm", "j.segm", "species",
                       "mhc.a", "mhc.b", "mhc.class", "antigen.epitope", "antigen.gene",
                       "antigen.species", "reference.id"),
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

#: Files compared as lines keyed on field 1, not by line position.
#:
#: The metadata files: one row per column of the table they describe, so the row identity is the
#: column name in field 1 and not the line number. Comparing them by position reported the single
#: inserted ``TCR_hash`` row as eight changed lines, because every row after it shifted down one.
#:
#: Their order is gated elsewhere, twice over: ``test_header_is_the_name_column_of_the_metadata``
#: asserts the name column is the data file's header, and :func:`_compare_table` fails a candidate
#: whose header is reordered. It is also why these are read as lines rather than parsed as a table:
#: the shipped ``web.method`` row has a space where a tab belongs, so the reference is ragged.
#:
#: ``latest-version.txt`` is here for the same reason and was compared by position until 2026-09-29.
#: It is a history with the newest release **prepended**, so every one of its 40 lines shifted down
#: and the comparison reported 40 changed lines and one added for a file that had gained exactly one
#: entry. It carries no tab, so field 1 is the entry itself, which is the key. The one
#: thing this cannot see - that line 1 is the *newest* release - is checked by `release.yml`'s
#: "Assert line 1 resolves" step and by `verify-latest.yml`, which compare it against the published
#: tag.
KEYED_LINES = ("vdjdb.meta.txt", "vdjdb.slim.meta.txt", "latest-version.txt")

#: Identity fields that live inside a JSON column rather than in a column of their own. They are
#: the seven ``meta.*`` members of :data:`CHUNK_DEDUP_KEY`: without them a study reporting one TCR
#: in several donors is one key group, and its rows pair arbitrarily.
JSON_KEYS: dict[str, tuple[str, tuple[str, ...]]] = {
    "vdjdb.txt": ("meta", ("study.id", "cell.subset", "subject.cohort", "subject.id",
                           "replica.id", "clone.id", "tissue")),
}

#: Surrogate keys: the values mean nothing on their own, only the partition they induce.
#:
#: ``complex.id`` is a counter allocated while walking the master table, so its values follow the
#: chunk order ``os.listdir`` happened to return. Comparing it by value made 185,868 of 284,546
#: rows "changed" in the first run -- 58 % of the table, none of it a difference in the data.
#: It is instead renumbered canonically on both sides before comparison, so a change in which
#: chains are grouped into one clone still shows up, and a reshuffle does not.
SURROGATE: dict[str, str] = {"vdjdb.txt": "complex.id", "vdjdb.slim.txt": "complex.id"}


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
    #: Whether a row-level comparison ran. When it did, its result is authoritative and a differing
    #: canonical digest is informational -- surrogate renumbering alone changes the digest.
    compared_rows: bool = False
    #: Declared rename id -> how many reference cells it rewrote in this file.
    renames: dict[str, int] = field(default_factory=dict)


@dataclass
class DiffReport:
    missing: list[str] = field(default_factory=list)
    added: list[str] = field(default_factory=list)
    files: list[FileReport] = field(default_factory=list)
    unattributed: list[CellDiff] = field(default_factory=list)
    rule_counts: dict[str, int] = field(default_factory=dict)
    rule_expected: dict[str, int] = field(default_factory=dict)
    row_deltas: dict[str, RowDelta] = field(default_factory=dict)
    rename_declared: list[str] = field(default_factory=list)
    #: ``file -> the instrument that gates it instead``. A member whose row-level comparison is not
    #: the instrument: it is still read, digested and counted, and the report says what does gate it,
    #: but its cells are not required to match a rule and its row buckets are not required to match a
    #: declaration. See :func:`load_elsewhere`.
    elsewhere: dict[str, str] = field(default_factory=dict)
    #: True when this report was loaded from a summary rather than computed, in which case every
    #: ``FileReport.cells`` is empty because the cells were not serialised. Anything that needs the
    #: per-cell detail must assert this is False rather than read an empty tuple as "no changes".
    cells_omitted: bool = False

    @property
    def rename_counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for f in self.files:
            for rid, n in f.renames.items():
                out[rid] = out.get(rid, 0) + n
        return out

    @property
    def stale_renames(self) -> list[str]:
        """Declared renames that matched nothing; a declaration outliving its data fails the run."""
        fired = self.rename_counts
        return sorted(r for r in self.rename_declared if not fired.get(r))

    @property
    def miscounted_rules(self) -> list[str]:
        """Rules that fired a different number of times than they declare."""
        return sorted(r for r, want in self.rule_expected.items()
                      if want >= 0 and self.rule_counts.get(r, 0) != want)

    @property
    def mismatched_row_deltas(self) -> list[str]:
        """Files whose measured row buckets differ from what ``[[row_delta]]`` declares.

        A named property rather than a branch inside :attr:`ok`, so :func:`render` can say which file
        and by how much. It could not before, and a run failing on this alone reported a FAIL with no
        finding anywhere in it - 0 unattributed cells, 0 miscounted rules, no stale rename - because
        the declared numbers were never printed. A new chunk is the common cause: its rows are rows
        the reference cannot contain, so all three totals move at once.
        """
        return sorted(f.name for f in self.files
                      if f.name not in self.elsewhere and self._delta_of(f) != self._measured_of(f))

    def _delta_of(self, f: FileReport) -> tuple[int, int]:
        d = self.row_deltas.get(f.name)
        return (d.removed, d.added) if d else (0, 0)

    @staticmethod
    def _measured_of(f: FileReport) -> tuple[int, int]:
        return (f.only_in_reference, f.only_in_candidate)

    @property
    def uncompared(self) -> list[str]:
        """Files judged on their canonical digest alone, and failing it.

        A row-level comparison not running is not itself a failure - a file with no identity key is
        compared by digest - but a digest that also differs leaves nothing said about the file, which
        is the one case where the comparison has no answer rather than a bad one.
        """
        return sorted(f.name for f in self.files
                      if f.name not in self.elsewhere and not f.compared_rows
                      and not f.canonical_equal)

    @property
    def ok(self) -> bool:
        return not (self.missing or self.added or self.unattributed or self.miscounted_rules
                    or self.stale_renames or self.mismatched_row_deltas or self.uncompared)


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

def _key_expr(name: str, df: pl.DataFrame, key: tuple[str, ...]) -> pl.Expr:
    parts: list[pl.Expr] = [pl.col(c) for c in key if c in df.columns]
    col, fields = JSON_KEYS.get(name, ("", ()))
    if col and col in df.columns:
        struct = pl.Struct([pl.Field(f, pl.Utf8) for f in fields])
        decoded = pl.col(col).str.json_decode(dtype=struct)
        parts += [decoded.struct.field(f).fill_null("") for f in fields]
    return pl.concat_str(parts, separator="\x1f")


def _contested(a: pl.DataFrame, b: pl.DataFrame) -> pl.Series:
    """Identity keys whose row multisets differ, found without materialising a single Python row.

    Materialising both tables costs 6.26 million tuples on ``vdjdb.txt`` and dominates the run.
    Instead, hash each row once in polars, count ``(key, row hash)`` pairs on both sides, and keep
    only the keys where some count disagrees. On a rebuild of the current corpus that is 14 keys of
    208,447, so the Python-level pairing below runs on ~0.007 % of the table.
    """
    def counts(df: pl.DataFrame) -> pl.DataFrame:
        return (df.group_by("__key", "__hash").len()
                  .rename({"len": "n"}))

    joined = counts(a).join(counts(b), on=["__key", "__hash"], how="full", coalesce=True)
    differing = joined.filter(
        pl.col("n").fill_null(0) != pl.col("n_right").fill_null(0)
    )
    return differing["__key"].unique()


def _rows_by_key(df: pl.DataFrame, keys: pl.Series) -> dict[str, list[tuple[str, ...]]]:
    """Materialise only the rows whose identity key is in ``keys``."""
    sub = df.filter(pl.col("__key").is_in(keys.implode()))
    out: dict[str, list[tuple[str, ...]]] = defaultdict(list)
    cols = [c for c in sub.columns if not c.startswith("__")]
    for k, row in zip(sub["__key"].to_list(), sub.select(cols).iter_rows(), strict=True):
        out[k].append(row)
    return out


def _canonical_surrogate(name: str, df: pl.DataFrame, col: str,
                         key: tuple[str, ...]) -> pl.DataFrame:
    """Renumber ``col`` so equal groupings get equal numbers, whatever order they were allocated in.

    Each group is labelled by the sorted row identities of its members; groups are then numbered by
    that label. ``0`` means "not grouped" (an unpaired chain) and is left alone.
    """
    label = (df.with_columns(_key_expr(name, df, key).alias("__k"))
               .with_row_index("__row")
               .group_by(col, maintain_order=True)
               .agg(pl.col("__k").sort().str.join("\x1e").alias("__label"),
                    pl.col("__row").min().alias("__i")))
    # `sort` must be stable and the tiebreak explicit: groups with an identical label are exactly
    # the ambiguous ones, and an unstable sort there made the comparison non-reproducible
    # (158 / 152 / 158 changed rows across three identical runs).
    order = (label.filter(pl.col(col) != "0")
                  .sort("__label", "__i", maintain_order=True)
                  .with_row_index("__n")
                  .with_columns((pl.col("__n") + 1).cast(pl.Utf8).alias("__new")))
    return (df.join(order.select(col, "__new"), on=col, how="left")
              .with_columns(pl.col("__new").fill_null("0").alias(col))
              .drop("__new"))


#: Fixed so a row hash is the same value in every run, on every host, and in every process.
#: polars seeds its hash from the value alone when given one, so a fixed seed makes the comparison
#: reproducible rather than merely repeatable.
_HASH_SEED = 20260925

#: Above this many unmatched rows in one key group, pair by sort order instead of by best match.
#: The quadratic search is affordable at 2-10 rows and not at 1,000.
_PAIR_LIMIT = 24


def _pair(left: list[tuple[str, ...]],
          right: list[tuple[str, ...]]) -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    """Pair rows that share a key but are not identical, minimising the reported difference.

    Sort order is the wrong correspondence: two rows of one study that differ in a field outside
    the identity key get crossed, and then every field that distinguishes them reports as changed.
    Greedy nearest-match pairs each row with the candidate it differs from least, so the report
    lists the smallest set of differences consistent with the data rather than an artifact of
    collation.
    """
    if not left or not right:
        return []
    if max(len(left), len(right)) > _PAIR_LIMIT:
        return list(zip(left, right, strict=False))
    free = list(range(len(right)))
    out = []
    for lrow in left:
        if not free:
            break
        j = min(free, key=lambda i: sum(a != b for a, b in zip(lrow, right[i], strict=True)))
        free.remove(j)
        out.append((lrow, right[j]))
    return out


def _compare_table(name: str, ref: bytes, cand: bytes,
                   renames: list[Rename] | None = None) -> FileReport:
    raw_r, can_r = _digests(ref)
    raw_c, can_c = _digests(cand)
    a, b = _read_table(ref), _read_table(cand)
    rename_counts: dict[str, int] = {}
    if renames:
        # The reference, never the candidate: the candidate is what we are measuring.
        a, rename_counts = _apply_renames(name, a, renames)

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

    surrogate = SURROGATE.get(name)
    if surrogate and surrogate in a.columns:
        a = _canonical_surrogate(name, a, surrogate, key)
        b = _canonical_surrogate(name, b, surrogate, key)

    cols = a.columns
    ident = _key_expr(name, a, key)
    a = a.with_columns(ident.alias("__key"),
                       pl.concat_str(cols, separator="\x1f", ignore_nulls=True).hash(seed=_HASH_SEED)
                       .alias("__hash"))
    b = b.with_columns(_key_expr(name, b, key).alias("__key"),
                       pl.concat_str(cols, separator="\x1f", ignore_nulls=True).hash(seed=_HASH_SEED)
                       .alias("__hash"))

    contested = _contested(a, b)
    ref_rows, cand_rows = _rows_by_key(a, contested), _rows_by_key(b, contested)
    only_ref = only_cand = changed = 0
    cells: list[CellDiff] = []

    for k in sorted(set(ref_rows) | set(cand_rows)):
        left, right = Counter(ref_rows.get(k, [])), Counter(cand_rows.get(k, []))
        common = left & right                      # identical rows, in multiset arithmetic
        lrows = sorted((left - common).elements())
        rrows = sorted((right - common).elements())
        for lrow, rrow in _pair(lrows, rrows):
            changed += 1
            cells.extend(CellDiff(name, c, lv, rv, k)
                         for c, lv, rv in zip(cols, lrow, rrow, strict=True) if lv != rv)
        only_ref += max(len(lrows) - len(rrows), 0)
        only_cand += max(len(rrows) - len(lrows), 0)

    return FileReport(name, raw_r == raw_c, can_r == can_c, a.height, b.height,
                      only_ref, only_cand, changed, tuple(cells), compared_rows=True,
                      renames=rename_counts)


def _compare_keyed_lines(name: str, ref: bytes, cand: bytes) -> FileReport:
    """Compare two metadata files by column name, reporting a whole differing line as one cell."""
    raw_r, can_r = _digests(ref)
    raw_c, can_c = _digests(cand)
    a = {line.split("\t")[0]: line for line in ref.decode().splitlines()}
    b = {line.split("\t")[0]: line for line in cand.decode().splitlines()}
    cells = tuple(CellDiff(name, "line", a[k], b[k], k)
                  for k in sorted(a.keys() & b.keys()) if a[k] != b[k])
    return FileReport(name, raw_r == raw_c, can_r == can_c, len(a), len(b),
                      len(a.keys() - b.keys()), len(b.keys() - a.keys()),
                      len(cells), cells, compared_rows=True)


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
    #: For a JSON column, the member that must be among the differing ones. Without it, one rule
    #: on `meta` covers every field inside it, so a regression in any of those fields is
    #: attributed to that rule.
    json_field: str | None = None

    def matches(self, c: CellDiff, json_fields: frozenset[str] | None = None) -> bool:
        """Does this rule cover ``c``?

        ``json_fields`` is :func:`_json_diff_fields` for ``c``, computed by the caller once per cell
        because every rule asks the same question of the same cell. Computing it here instead
        parsed the two JSON cells once per *rule*: 11,390,336 calls for 201,481 cells, which an
        `lru_cache(maxsize=4096)` on the parse was covering for rather than fixing.
        """
        if not ((self.file is None or self.file == c.file)
                and (self.column is None or self.column == c.column)
                and (self.from_ is None or self.from_ == c.old)
                and (self.to is None or self.to == c.new)):
            return False
        if self.json_field is None:
            return True
        return self.json_field in (json_fields if json_fields is not None
                                   else _json_diff_fields(c))


@dataclass(frozen=True, slots=True)
class Rename:
    """A declared value rewrite, applied to the reference before rows are keyed.

    A nomenclature correction to a column that is part of the identity key produces no changed cell:
    it removes a row on one side and adds one on the other, and the cell-level machinery has nothing
    to attribute. Rewriting the reference first restores the match, so the comparison still measures
    what else moved.

    The rewrite itself is not verified here; it is reviewed as a diff of ``rules/`` and counted by
    the build's own `nomenclature.tsv` report. A rename that matches nothing fails the run, so a
    declaration cannot outlive the data it describes.
    """

    columns: tuple[str, ...]
    from_: str
    to: str
    files: tuple[str, ...] = ()          # empty = every table
    #: The organism the rewrite is for. Empty = every organism, which is safe only while ``from_``
    #: is a spelling nobody uses correctly, and stops being safe the moment it is a real gene
    #: somewhere else. ``TRAV14-1*01`` is one human record's respelling of ``TRAV14/DV4`` and the
    #: correct IMGT name of 79 mouse records in the same release: unscoped it rewrote all 80 and
    #: keyed 79 mouse rows against a gene that is not theirs (#671).
    species: str = ""
    #: Optional evidence: rewrite only rows where one of these columns contains ``when_contains``.
    #: Without it a rename must be injective on value alone, which an evidence-based correction is
    #: not -- the 1,047 TRAJ24 records the CDR3 identifies as ``*02`` ship the same ``TRAJ24*01`` as
    #: the 38 it does not, so an unconditional rename would rewrite both.
    when_columns: tuple[str, ...] = ()
    when_contains: str = ""
    #: Exact alternative to ``when_contains``. Required when the evidence column holds peptides: a
    #: 9-mer epitope is a substring of a 10-mer one (`SPRWYFYYL` inside `LSPRWYFYYL`), so a
    #: substring predicate would fire on the wrong rows.
    when_equals: str = ""

    @property
    def id(self) -> str:
        base = f"{self.from_} -> {self.to}"
        evidence = self.when_contains or self.when_equals
        return f"{base} [{evidence}]" if evidence else base

    @property
    def conditional(self) -> bool:
        """True when the rewrite needs a predicate, so it cannot go in a plain value mapping."""
        return bool(self.species or (self.when_columns
                                     and (self.when_contains or self.when_equals)))

    def predicate(self, snap: dict[str, str]) -> pl.Expr:
        """The full condition, read from ``snap``: the species scope and the evidence, if any."""
        parts = []
        if self.species and "species" in snap:
            parts.append(pl.col(snap["species"]) == self.species)
        evidence = [snap[c] for c in self.when_columns if c in snap]
        if evidence and (self.when_contains or self.when_equals):
            parts.append(self.evidence_expr(evidence))
        return pl.all_horizontal(*parts) if parts else pl.lit(True)

    def evidence_expr(self, columns: list[str]) -> pl.Expr:
        """True where one of ``columns`` contains the declared evidence."""
        if self.when_equals:
            return pl.any_horizontal(*[pl.col(c) == self.when_equals for c in columns])
        return pl.any_horizontal(
            *[pl.col(c).str.contains(self.when_contains, literal=True) for c in columns])

    def applies_to(self, file: str) -> bool:
        return not self.files or file in self.files


def load_renames(path: Path) -> list[Rename]:
    if not path.exists():
        return []
    raw = tomllib.loads(path.read_text())
    return [Rename(columns=tuple(c.strip() for c in r["columns"].split(",")),
                   from_=r["from"], to=r["to"],
                   files=tuple(f.strip() for f in r.get("files", "").split(",") if f.strip()),
                   species=r.get("species", ""),
                   when_columns=tuple(c.strip() for c in r.get("when_columns", "").split(",")
                                      if c.strip()),
                   when_contains=r.get("when_contains", ""),
                   when_equals=r.get("when_equals", ""))
            for r in raw.get("rename", [])]


def _apply_renames(name: str, df: pl.DataFrame,
                   renames: list[Rename]) -> tuple[pl.DataFrame, dict[str, int]]:
    """Rewrite the reference's declared values, and count what each declaration matched.

    Simultaneously, per column, in one pass. Applying them one after another chains them: with
    ``A -> B`` and ``B -> C`` declared, a cell that was already ``B`` comes out ``C``. Chaining
    moved mouse ``TRAV6-1*01`` rows onto ``TRAV6-7/DV9*01`` with no error and turned a clean
    comparison into 6,334 phantom unmatched rows, in the reference, where nothing had changed.
    """
    counts: dict[str, int] = {}
    per_column: dict[str, dict[str, str]] = {}
    conditional: list[Rename] = []
    masks: list[tuple[str, pl.Expr]] = []
    for r in renames:
        if not r.applies_to(name):
            continue
        cols = [c for c in r.columns if c in df.columns]
        if not cols:
            continue
        mask = pl.any_horizontal(*[pl.col(c) == r.from_ for c in cols])
        if r.conditional:
            if r.species and "species" not in df.columns:
                continue
            evidence = [c for c in r.when_columns if c in df.columns]
            if r.when_columns and not evidence:
                continue
            mask = mask & r.predicate({c: c for c in df.columns})
            conditional.append(r)
        else:
            for c in cols:
                per_column.setdefault(c, {})[r.from_] = r.to
        masks.append((r.id, mask))

    # **One `select` for every declared count, not one per rename.** `df.select(mask.sum()).item()`
    # inside the loop above was 353 independent query plans over the full table, and with the
    # per-column `with_columns` below it made this function 875 `collect()` calls for 32.4 s of a
    # 40.1 s `vdjdb diff`. Batched: 8.96 s, and `summary_json` is byte-identical.
    if masks:
        totals = df.select(*[m.sum().alias(f"__n{i}") for i, (_, m) in enumerate(masks)]).row(0)
        for (rid, _), n in zip(masks, totals, strict=True):
            counts[rid] = counts.get(rid, 0) + int(n or 0)

    if per_column:
        df = df.with_columns(*[pl.col(c).replace(m) for c, m in sorted(per_column.items())])
    if conditional:
        # Every conditional is evaluated against a snapshot taken before any of them apply, for the
        # same reason the unconditional ones share one mapping: otherwise they chain.
        # Both the target and the evidence are read from the snapshot, so a pair of renames that
        # *swap* two columns works: without it the second would test a column the first has already
        # rewritten, and the order of declaration would decide the answer.
        touched = {c for r in conditional
                   for c in (*r.columns, *r.when_columns, "species") if c in df.columns}
        snap = {c: f"__snap\x1f{c}" for c in sorted(touched)}
        df = df.with_columns(*[pl.col(c).alias(s) for c, s in snap.items()])
        # The column set is a **set**: a column named by several conditional renames was rewritten
        # once per rename, to the same value each time, paying a full materialisation for each. And
        # the expressions are applied in one `with_columns`, not one per column.
        rewritten = []
        for c in sorted({c for r in conditional for c in r.columns if c in snap}):
            expr = pl.col(c)
            for r in conditional:
                if c not in r.columns:
                    continue
                expr = pl.when((pl.col(snap[c]) == r.from_) & r.predicate(snap)
                               ).then(pl.lit(r.to)).otherwise(expr)
            rewritten.append(expr.alias(c))
        df = df.with_columns(*rewritten).drop(list(snap.values()))
    return df, counts


def _json_members(text: str) -> tuple[tuple[str, str], ...]:
    try:
        return tuple((k, repr(v)) for k, v in json.loads(text).items())
    except (ValueError, AttributeError):
        return ()


def _json_diff_fields(c: CellDiff) -> frozenset[str]:
    """The member names that differ between two JSON cells."""
    a, b = dict(_json_members(c.old)), dict(_json_members(c.new))
    if not a and not b:
        return frozenset()
    return frozenset(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))


@dataclass(frozen=True, slots=True)
class RowDelta:
    """A declared, measured change in which rows a file contains.

    A nomenclature correction to a column that is part of a file's grouping key removes one row and
    adds another -- there is no cell to attribute. Declaring the counts records that without
    weakening the gate: an undeclared row delta still fails.
    """

    file: str
    added: int
    removed: int
    note: str = ""


def load_row_deltas(path: Path) -> dict[str, RowDelta]:
    """One declaration per file, and a second one for the same file is an error rather than a winner.

    The dict comprehension this replaces let a later ``[[row_delta]]`` silently overwrite an earlier
    one for the same file. That is the worst possible failure for this file: a curator declaring the
    rows their chunk adds would erase the declaration describing the code deviations, and the
    comparison would still read PASS on an accounting that no longer covers both. Raising makes the
    collision a build failure with both line numbers' worth of context in the message.
    """
    if not path.exists():
        return {}
    raw = tomllib.loads(path.read_text())
    out: dict[str, RowDelta] = {}
    for d in raw.get("row_delta", []):
        name = d["file"]
        if name in out:
            raise ValueError(
                f"{path}: two [[row_delta]] declarations for {name!r}. There is one per file, so "
                f"combine them into a single added/removed pair and say in `note` what each "
                f"component is. Existing: added={out[name].added} removed={out[name].removed}; "
                f"second: added={int(d.get('added', 0))} removed={int(d.get('removed', 0))}")
        out[name] = RowDelta(name, int(d.get("added", 0)), int(d.get("removed", 0)),
                             d.get("note", ""))
    return out


def load_members(path: Path) -> tuple[frozenset[str], frozenset[str]]:
    """``(declared added, declared removed)`` bundle members, from ``[members]``.

    A release that ships a file the reference does not contain is a declared change of shape, not an
    unexplained one - the TCREMP motif tables are the first of them - and until this existed the
    file-set pass had no way to say so: any new member failed the comparison outright, so the only
    way to run it was to name the reference's members with ``--only`` and never compare the rest.
    An **undeclared** new or missing member still fails, which is the half that matters.
    """
    if not path.exists():
        return frozenset(), frozenset()
    m = tomllib.loads(path.read_text()).get("members", {})
    return frozenset(m.get("added", ())), frozenset(m.get("removed", ()))


def load_elsewhere(path: Path) -> dict[str, str]:
    """``[measured_elsewhere]``: member -> the instrument that gates it instead of this comparison.

    Three members of the legacy bundle cannot be judged by keying their rows against the reference,
    and saying so is better than the alternative that was in place, which was to name the five
    comparable members with ``--only`` in a workflow file and never compare the other seven at all.
    Measured on the assembled 2026-09-29 legacy zip: that left **101,877 unattributed cells**, every
    one of them in a motif file or in ``latest-version.txt``, and nothing in the repository said
    whether that was expected.

    A file named here is still read, digested, row-counted and printed. What it is exempt from is the
    requirement that every changed cell match a rule and every row bucket match a declaration - and
    the value says which instrument takes that job, so the exemption names its replacement rather
    than removing a check.
    """
    if not path.exists():
        return {}
    raw = tomllib.loads(path.read_text()).get("measured_elsewhere", {})
    return {str(k): str(v) for k, v in raw.items()}


def load_rules(path: Path) -> list[Rule]:
    if not path.exists():
        return []
    raw = tomllib.loads(path.read_text())
    return [Rule(id=r["id"], file=r.get("file"), column=r.get("column"),
                 from_=r.get("from"), to=r.get("to"), rows=int(r.get("rows", -1)),
                 json_field=r.get("json_field"))
            for r in raw.get("rule", [])]


# ---------------------------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------------------------

def summary_json(report: DiffReport) -> str:
    """The report without its per-cell detail, as JSON. Small enough to be an artifact.

    **The cells are deliberately left out.** ``dataclasses.asdict`` on a real report is **211 MB**,
    because the declared differences are 485,130 changed cells across five files, and serialising then
    parsing that costs more than recomputing the comparison. Every scalar the gate and the release
    tests read is here; ``unattributed`` is kept in full because it is the failure evidence and is
    empty on a PASS by construction.

    ``cells_omitted`` is set on load so nothing can read an empty ``cells`` tuple as "no changes"
    (``ROADMAP_local.md`` section 57.3).
    """
    import dataclasses
    import json

    d = dataclasses.asdict(report)
    for f in d["files"]:
        f["cells"] = []
    d["cells_omitted"] = True
    return json.dumps(d, indent=1, sort_keys=True)


def summary_from_json(text: str) -> DiffReport:
    """Rebuild a :class:`DiffReport` from :func:`summary_json`. ``cells`` is empty by construction."""
    import json

    d = json.loads(text)
    return DiffReport(
        missing=d["missing"], added=d["added"],
        files=[FileReport(**{**f, "cells": (), "renames": f.get("renames", {})})
               for f in d["files"]],
        unattributed=[CellDiff(**c) for c in d["unattributed"]],
        rule_counts=d["rule_counts"], rule_expected=d["rule_expected"],
        row_deltas={k: RowDelta(**v) for k, v in d["row_deltas"].items()},
        rename_declared=d["rename_declared"],
        # Without this the round trip drops the exemptions and `ok` reads False on a report that
        # passed: the release tests read the build's JSON rather than recomputing the comparison, so
        # a field missing here is a field the gate does not have.
        elsewhere=d.get("elsewhere", {}), cells_omitted=True)


def diff(reference: Path, candidate: Path, rules_path: Path | None = None,
         only: Iterable[str] | None = None) -> DiffReport:
    """Compare two bundles and attribute every changed cell to a declared rule.

    ``only`` restricts the comparison to the named members, for candidates that are deliberately
    partial -- the assembly stage produces three tables, and the motif and dashboard members arrive
    from later stages.
    """
    ref, cand = Bundle(reference), Bundle(candidate)
    wanted = set(only) if only else None
    ref_names = [n for n in ref.names if wanted is None or n in wanted]
    cand_names = [n for n in cand.names if wanted is None or n in wanted]
    rules = load_rules(rules_path) if rules_path else []
    renames = load_renames(rules_path) if rules_path else []
    declared_added, declared_removed = load_members(rules_path) if rules_path else (frozenset(),
                                                                                    frozenset())
    report = DiffReport(
        missing=[n for n in ref_names if n not in cand_names and n not in declared_removed],
        added=[n for n in cand_names if n not in ref_names and n not in declared_added],
        rule_expected={r.id: r.rows for r in rules},
        row_deltas=load_row_deltas(rules_path) if rules_path else {},
        rename_declared=sorted({r.id for r in renames}),
        elsewhere=load_elsewhere(rules_path) if rules_path else {},
    )

    for name in ref_names:
        # `cand_names`, not `report.missing`: a member declared under `[members] removed` is absent
        # from the candidate and absent from the missing list, and reading it raised a KeyError.
        if name not in cand_names:
            continue
        a, b = ref.read_bytes(name), cand.read_bytes(name)
        if name in KEYS or name in WIDTHS:
            report.files.append(_compare_table(name, a, b, renames))
        elif name in KEYED_LINES:
            report.files.append(_compare_keyed_lines(name, a, b))
        else:
            report.files.append(_compare_opaque(name, a, b))

    for f in report.files:
        if f.name in report.elsewhere:
            continue
        # Any rule carrying a `json_field` asks which members of this cell differ, and every rule
        # asks it of the same cell, so it is computed once here. Only cells a JSON rule could match
        # pay for it at all.
        json_columns = {r.column for r in rules if r.json_field is not None}
        for cell in f.cells:
            fields = _json_diff_fields(cell) if cell.column in json_columns else None
            for rule in rules:
                if rule.matches(cell, fields):
                    report.rule_counts[rule.id] = report.rule_counts.get(rule.id, 0) + 1
                    break
            else:
                report.unattributed.append(cell)
    return report


def render(report: DiffReport) -> str:
    """The release-notes table, which ``--report`` writes."""
    out = ["# Comparison against the reference release", ""]
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

    if report.elsewhere:
        out += ['', '## Gated by another instrument', '',
                'These members are read, digested and counted above, and their cells are **not**',
                'required to match a rule here. Keying their rows against the reference answers no',
                'question a reader has: a `cid` carries a position in a sorted list, so one',
                'renumbered cluster reads as the entire file replaced.', '',
                '| Member | What gates it |', '|---|---|']
        for name, how in sorted(report.elsewhere.items()):
            out.append(f'| `{name}` | {how} |')

    if report.rename_declared:
        out += ['', '## Declared renames', '',
                'Applied to the reference before keying, so a nomenclature correction does not read',
                'as a lost row and a found one.', '',
                '| Rename | Reference cells rewritten |', '|---|---|']
        fired = report.rename_counts
        for rid in report.rename_declared:
            out.append(f'| `{rid}` | {fired.get(rid, 0)} |')
        if report.stale_renames:
            out += ['', '**Stale renames (matched nothing):** '
                    + ', '.join(f'`{r}`' for r in report.stale_renames)]

    if report.rule_expected:
        out += ["", "## Declared rules", "", "| Rule | Declared | Fired |", "|---|---|---|"]
        for rid, want in sorted(report.rule_expected.items()):
            got = report.rule_counts.get(rid, 0)
            flag = "" if want < 0 or want == got else "  **MISCOUNT**"
            out.append(f"| `{rid}` | {'--' if want < 0 else want} | {got}{flag} |")

    if report.row_deltas or report.mismatched_row_deltas:
        out += ["", "## Declared row buckets", "",
                "A row whose identity changed leaves the reference bucket and arrives in the candidate",
                "one, so these two counts are declared per file the way a rule's count is.", "",
                "| File | Declared only ref | Measured | Declared only cand | Measured |",
                "|---|---|---|---|---|"]
        for f in sorted(report.files, key=lambda f: f.name):
            want_ref, want_cand = report._delta_of(f)
            if not (report.row_deltas.get(f.name) or f.only_in_reference or f.only_in_candidate):
                continue
            # A member gated by another instrument is counted here and not flagged: bolding a bucket
            # this comparison is not judging reads as a failure it then does not fail on.
            gated = f.name not in report.elsewhere
            def flag(want: int, got: int, gated: bool = gated) -> str:
                return f"{got}" if not gated or want == got else f"**{got}**"
            want = (f"{want_ref}", f"{want_cand}") if gated else ("-", "-")
            out.append(f"| `{f.name}` | {want[0]} | {flag(want_ref, f.only_in_reference)} | "
                       f"{want[1]} | {flag(want_cand, f.only_in_candidate)} |")
        if report.mismatched_row_deltas:
            out += ["", "**MISMATCH** in bold above: "
                    + ", ".join(f"`{n}`" for n in report.mismatched_row_deltas)
                    + ". Re-measure and update the `[[row_delta]]` block in "
                      "`rules/expected_diffs.toml`, extending its `note` with the reason - a landing "
                      "chunk names the chunk and its record count."]

    if report.uncompared:
        out += ["", "**Compared by digest only, and the digest differs:** "
                + ", ".join(f"`{n}`" for n in report.uncompared)
                + ". No row-level statement was made about these files."]

    if report.unattributed:
        by_col = Counter((c.file, c.column) for c in report.unattributed)
        out += ["", f"## Unattributed differences: {len(report.unattributed)} cells", "",
                "| File | Column | Cells | Example |", "|---|---|---|---|"]
        for (fname, col), n in by_col.most_common(40):
            ex = next(c for c in report.unattributed if c.file == fname and c.column == col)
            out.append(f"| `{fname}` | `{col}` | {n} | `{ex.old[:40]}` → `{ex.new[:40]}` |")

    out += ["", f"**Verdict: {'PASS' if report.ok else 'FAIL'}**"]
    return "\n".join(out) + "\n"
