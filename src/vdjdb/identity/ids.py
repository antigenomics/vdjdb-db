"""Stable VDJdb record identifiers.

A VDJdb record has no identifier today. It is a row in a chunk file, and the only thing that names it
is its own content. So when a curator fixes a typo in a CDR3, the old record disappears and a new one
appears, and nothing downstream can tell that from a deletion plus an addition. Citations to
vdjdb.com, structure links keyed on ``TCR_hash``, and evidence accumulated across releases all break,
with no error.

Identity here is therefore two-level:

``record_id``
    Opaque, assigned once, never reused, and stable across content changes. This is what external
    references point at. Format ``VDJDB<10 digits>``, allocated monotonically.

``content_hash``
    sha256 over the record's canonical content. Changes whenever anything about the record changes.
    This is what detects a change; it is never an identifier.

Plus a lifecycle (``active`` / ``amended`` / ``retired``) and provenance back to the chunk file, the
row within it, and the commit that last touched it.

Matching cascade. On each build, records from ``chunks/`` are reconciled against the committed
registry in three passes, most specific first:

1. **Exact** natural key -> the same ``record_id``, unchanged.
2. **Amendment**: within the same chunk and the same ``reference.id``, a record whose natural key
   differs from an unmatched registry entry in exactly one field is treated as that entry,
   amended. This makes a typo fix traceable instead of a delete-plus-insert. Ambiguous
   matches (two equally good candidates) are refused and fall through -- a wrong link is worse than
   a new id.
3. **New**: anything still unmatched gets a fresh id. Registry entries still unmatched afterwards
   are ``retired``, carrying the release in which they were last seen.

The registry is a committed TSV sorted by ``record_id`` so that a curation PR shows exactly which
records were added, amended or retired, as a reviewable diff.
"""
from __future__ import annotations

import gzip
import hashlib
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import IO, cast

import polars as pl

from ..schema import ALL_COLUMNS, CHUNK_DEDUP_KEY, LEGACY_CHUNK_DEDUP_KEY

#: The fields that identify a record. A change in any of these is an amendment; a change anywhere
#: else is re-annotation that keeps the same record.
#:
#: ``CHUNK_DEDUP_KEY`` plus the chunk. A chunk is one paper, so two rows in two chunks are
#: independent reports rather than one record seen twice, however identical their fields (CLAUDE.md,
#: the data model). Deduplication is within a chunk; identity is per curated line.
#:
#: Two earlier versions were wrong in opposite directions. A narrower key that stopped at
#: ``reference.id`` collided on 20,769 of 192,753 records -- one paper reporting the same TCR
#: against the same epitope in several donors is several records. Dropping ``chunk.file`` merged
#: 19 pairs that are two papers' independent reports, the signal phase 11 tunes against.
NATURAL_KEY: tuple[str, ...] = (*CHUNK_DEDUP_KEY, "chunk.file")
LEGACY_NATURAL_KEY: tuple[str, ...] = (*LEGACY_CHUNK_DEDUP_KEY, "chunk.file")

ID_PREFIX = "VDJDB"
ID_DIGITS = 10

REGISTRY_COLUMNS: tuple[str, ...] = (
    "record_id",
    "state",
    "natural_key_hash",
    "content_hash",
    "chunk_file",
    "chunk_row",
    "first_seen_release",
    "first_seen_commit",
    "last_seen_release",
    "last_modified_release",
    "last_modified_commit",
    "amendment_count",
    "amended_from_key_hash",
    "replaced_by",
    "note",
)

#: The receptor a natural key names, by position in :data:`NATURAL_KEY`. Two keys agreeing on all
#: six describe the same TCR, whatever moved around it -- the test :func:`_successor` applies before
#: it will link a retired id to a new one.
RECEPTOR_FIELDS: tuple[int, ...] = tuple(
    NATURAL_KEY.index(c) for c in
    ("cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta"))


class RecordState(StrEnum):
    ACTIVE = "active"
    RETIRED = "retired"


def format_id(n: int) -> str:
    return f"{ID_PREFIX}{n:0{ID_DIGITS}d}"


def parse_id(rid: str) -> int:
    if not rid.startswith(ID_PREFIX):
        raise ValueError(f"not a VDJdb record id: {rid!r}")
    return int(rid[len(ID_PREFIX):])


def _hash(parts: list[str]) -> str:
    """sha256 over tab-joined parts. Units separator guards against a field containing a tab."""
    return hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()


def natural_key(row: dict[str, object]) -> str:
    """Hash of the identifying fields. Two records with the same value are the same record."""
    return _hash([str(row.get(c) or "").strip() for c in NATURAL_KEY])


def canonical_content_hash(row: dict[str, object]) -> str:
    """Hash of every build-relevant field, in a fixed order.

    Changes whenever anything about the record changes, including method and meta.
    The natural key answers "is this the same record"; the content hash
    answers "has it changed".
    """
    return _hash([str(row.get(c) or "").strip() for c in ALL_COLUMNS])


@dataclass(slots=True)
class _Entry:
    record_id: str
    state: str
    natural_key_hash: str
    content_hash: str
    chunk_file: str
    chunk_row: int
    first_seen_release: str
    first_seen_commit: str
    last_seen_release: str
    last_modified_release: str
    last_modified_commit: str
    amendment_count: int
    amended_from_key_hash: str
    #: The id that took over when this one retired because two key fields moved at once, so a
    #: consumer holding the old id can follow it forward (`ROADMAP.md` §10.4, #693). Empty on an
    #: active record, and empty on a retirement that is a genuine deletion.
    replaced_by: str
    note: str


@dataclass
class ReconcileReport:
    """What happened in one reconciliation. This is the reviewable output."""

    unchanged: int = 0
    annotated: int = 0          # same natural key, different content
    amended: list[tuple[str, str, str, str]] = field(default_factory=list)  # id, field, old, new
    added: list[str] = field(default_factory=list)
    retired: list[str] = field(default_factory=list)
    #: (record_id, chunk_file) for rows that are a second submission of an existing record.
    replicated: list[tuple[str, str]] = field(default_factory=list)

    def summary(self) -> str:
        return (
            f"unchanged {self.unchanged}  annotated {self.annotated}  "
            f"amended {len(self.amended)}  added {len(self.added)}  retired {len(self.retired)}  "
            f"replicated {len(self.replicated)}"
        )


class IdentityRegistry:
    """The committed record registry."""

    def __init__(self, entries: list[_Entry] | None = None) -> None:
        self._by_key: dict[str, _Entry] = {}
        self._next = 1
        for e in entries or []:
            self._by_key[e.natural_key_hash] = e
            self._next = max(self._next, parse_id(e.record_id) + 1)

    # ---- persistence -------------------------------------------------------

    @classmethod
    def load(cls, path: Path) -> IdentityRegistry:
        if not path.exists():
            return cls()
        # Empty string is the only missing marker (CLAUDE.md), and the polars kwarg for that has
        # been renamed across versions -- normalise after the read instead.
        df = pl.read_csv(path, separator="\t", infer_schema=False, quote_char=None).fill_null("")
        missing = set(REGISTRY_COLUMNS) - set(df.columns) - {"replaced_by"}
        if missing:
            raise ValueError(f"{path} is missing registry columns: {sorted(missing)}")
        # Older inputs predate the successor column. Preserve their empty successor values.
        if "replaced_by" not in df.columns:
            df = df.with_columns(pl.lit("").alias("replaced_by"))
        df = df.select(REGISTRY_COLUMNS).with_columns(
            pl.col("chunk_row", "amendment_count").str.strip_chars().replace("", "0").cast(pl.Int64)
        )
        entries = [_Entry(*row) for row in df.iter_rows()]
        return cls(entries)

    def save(self, path: Path) -> None:
        rows = sorted(self._by_key.values(), key=lambda e: parse_id(e.record_id))
        df = pl.DataFrame(
            {c: [str(getattr(e, c)) for e in rows] for c in REGISTRY_COLUMNS}
        ) if rows else pl.DataFrame({c: pl.Series(c, [], pl.String) for c in REGISTRY_COLUMNS})
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            # Fixed metadata makes the committed input identical across updates and filenames.
            with (path.open("wb") as raw,
                  gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed):
                df.write_csv(cast(IO[bytes], compressed), separator="\t", quote_style="never",
                             line_terminator="\n")
        else:
            df.write_csv(path, separator="\t", quote_style="never", line_terminator="\n")

    # ---- lookup ------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._by_key)

    def active(self) -> list[_Entry]:
        return [e for e in self._by_key.values() if e.state == RecordState.ACTIVE]

    def get(self, key_hash: str) -> _Entry | None:
        return self._by_key.get(key_hash)

    def _allocate(self) -> str:
        rid = format_id(self._next)
        self._next += 1
        return rid

    def to_frame(self) -> pl.DataFrame:
        rows = sorted(self._by_key.values(), key=lambda e: parse_id(e.record_id))
        return pl.DataFrame({c: [str(getattr(e, c)) for e in rows] for c in REGISTRY_COLUMNS})


def _key_fields(row: dict[str, object]) -> tuple[str, ...]:
    return tuple(str(row.get(c) or "").strip() for c in NATURAL_KEY)


def _single_field_difference(a: tuple[str, ...], b: tuple[str, ...]) -> int | None:
    """Index of the one differing field, or None if zero or more than one differ."""
    diff = [i for i, (x, y) in enumerate(zip(a, b, strict=True)) if x != y]
    return diff[0] if len(diff) == 1 else None


def reconcile(
    records: pl.DataFrame,
    registry: IdentityRegistry,
    *,
    release: str,
    commits: dict[str, str] | None = None,
) -> tuple[pl.DataFrame, IdentityRegistry, ReconcileReport]:
    """Assign a ``record_id`` to every row of ``records``, updating ``registry`` in place.

    ``records`` must have ``ALL_COLUMNS`` plus ``chunk.file`` and ``chunk.row``.
    ``commits`` maps a chunk filename to the commit that last touched it.
    """
    commits = commits or {}
    report = ReconcileReport()

    rows = records.to_dicts()
    keys = [_key_fields(r) for r in rows]
    key_hashes = [natural_key(r) for r in rows]
    content_hashes = [canonical_content_hash(r) for r in rows]

    # Migrate the old key only when both content and the original source row agree.
    # Newly retained metadata variants must not take the ID of the earlier retained row.
    for row, kh, content in zip(rows, key_hashes, content_hashes, strict=True):
        old = _hash([str(row.get(c) or "").strip() for c in LEGACY_NATURAL_KEY])
        entry = registry.get(old)
        if (entry is not None and entry.content_hash == content
                and entry.chunk_row == int(row.get("chunk.row") or 0)):
            del registry._by_key[old]
            entry.natural_key_hash = kh
            registry._by_key[kh] = entry

    matched_entries: set[str] = set()
    assigned: list[str | None] = [None] * len(rows)

    # Pass 1 -- exact natural key. Several rows may share one key: the same record submitted in two
    # chunk is one record; ``chunk.file`` is part of the key, so two chunks never collide.
    seen_key: dict[str, int] = {}
    for i, kh in enumerate(key_hashes):
        e = registry.get(kh)
        if e is not None:
            first = seen_key.get(kh)
            if first is not None:
                assigned[i] = assigned[first]
                report.replicated.append((str(assigned[first]), str(rows[i].get("chunk.file") or "")))
                continue
            seen_key[kh] = i
            assigned[i] = e.record_id
            matched_entries.add(e.record_id)
            e.state = RecordState.ACTIVE
            e.last_seen_release = release
            e.chunk_file = str(rows[i].get("chunk.file") or e.chunk_file)
            e.chunk_row = int(rows[i].get("chunk.row") or 0)
            if e.content_hash != content_hashes[i]:
                e.content_hash = content_hashes[i]
                e.last_modified_release = release
                e.last_modified_commit = commits.get(e.chunk_file, "")
                report.annotated += 1
            else:
                report.unchanged += 1

    # Pass 2 -- amendment: same chunk and reference, exactly one differing key field.
    unmatched_rows = [i for i, a in enumerate(assigned) if a is None]
    if unmatched_rows:
        leftovers = [e for e in registry.active() if e.record_id not in matched_entries]
        # Bucket candidates by chunk_file, so the scan stays local.
        #
        # ⚠ Not by `(chunk_file, reference.id)`, which is what this did until #685. The bucket was
        # built from the *previous* build's reference and read with the *new* one, so when
        # `reference.id` was the field that moved the two never agreed, the bucket came back empty,
        # and the single-field rule this pass exists for was never reached: the record was retired
        # and a new id minted, silently. Closing the space in `PMID: 34433824` retired 22 published
        # ids that way.
        #
        # The reference was not keeping the scan local either. A chunk is one paper, so 213 of the
        # 231 chunks carry exactly one `reference.id`; the 18 that carry more are small, the two
        # largest being `small_datasets_2026-05-29.tsv` (1,044 rows, 237 references) and
        # `PDB_Database.tsv` (370 / 209), against a 29,715-row single-reference chunk the second
        # component did nothing for. The chunk does the localising.
        by_bucket: dict[str, list[_Entry]] = defaultdict(list)
        by_position: dict[tuple[str, int], list[_Entry]] = defaultdict(list)
        for e in leftovers:
            if _unpack_note(e.note) is not None:
                by_bucket[e.chunk_file].append(e)
                by_position[e.chunk_file, e.chunk_row].append(e)

        # The registry does not store the raw key, only its hash, so amendment matching needs the
        # previous build's key fields. They are recorded in `note` as the reference id plus the
        # packed key; see `_pack_note`. Entries written by older versions do not match.
        # Every row number each chunk still has, for the "no line moved" test below. Built from all
        # rows of the chunk, not only the unmatched ones: what makes a row number an identity is that
        # the file still has a row there, whether or not that row needs matching.
        rows_present: dict[str, set[int]] = defaultdict(set)
        for r in rows:
            rows_present[str(r.get("chunk.file") or "")].add(int(r.get("chunk.row") or 0))

        intact_chunks = {name for name, entries in by_bucket.items()
                         if all(e.chunk_row in rows_present[name] for e in entries)}
        for i in unmatched_rows:
            row, key = rows[i], keys[i]
            chunk = str(row.get("chunk.file") or "")
            bucket = by_bucket.get(chunk, [])
            candidates: list[tuple[_Entry, int]] = []
            # When every old position remains, the existing ambiguity rule below selects the
            # same-row candidate. Resolve that case directly instead of scanning a large chunk.
            if chunk in intact_chunks:
                for e in by_position.get((chunk, int(row.get("chunk.row") or 0)), []):
                    if e.record_id not in matched_entries:
                        prev = _unpack_note(e.note)
                        assert prev is not None
                        d = _single_field_difference(key, prev)
                        if d is not None:
                            candidates.append((e, d))
            if len(candidates) == 1:
                bucket = []
            else:
                candidates = []
            for e in bucket:
                if e.record_id in matched_entries:
                    continue
                prev = _unpack_note(e.note)
                if prev is None:
                    continue
                d = _single_field_difference(key, prev)
                if d is not None:
                    candidates.append((e, d))
            if len(candidates) > 1:
                # Two records of one chunk can both be one field from this row, and the field's value
                # is then not enough to say which. `chunk.row` is - but only when no line moved, because
                # a deleted line shifts every number after it and the number alone would be a guess.
                #
                # What separates the two is whether the chunk still has a row at every number the
                # candidates sit on. An edit in place leaves all of them there, so matching by number
                # is an identity: the same line of the same file. A removed line takes its number with
                # it, the test fails, the ambiguity stands, and pass 3 mints a new id - which is what
                # `test_ambiguous_amendment_is_refused` asks for and what
                # `test_the_row_tiebreak_is_refused_when_a_line_moved` asserts from the other side.
                #
                # Two cases needed it, and the first is why the condition cannot be narrower. Repairing
                # the 4,324 junctions of vdjdb-db#646 moves that many keys at once, so a chunk's bucket
                # holds hundreds of unmatched rows and each row has its own two or three candidates: a
                # bijection between the bucket's rows and *one row's* candidates never holds, and the
                # first version of this test therefore refused every one of them - 2,066 published ids
                # retired and re-minted. The second is `menon_etal_2024.tsv` rows 26 and 27, whose
                # `TRBV5-3;TRBV5-5;TRBV5-8` and `TRBV5-3;TRBV5-8` both moved one field when `;` became
                # `,` (`52cb4e2`), costing `VDJDB0000187889` and `...890` their ids.
                present = rows_present.get(str(row.get("chunk.file") or ""), set())
                if all(c[0].chunk_row in present for c in candidates):
                    candidates = [c for c in candidates
                                  if c[0].chunk_row == int(row.get("chunk.row") or 0)]
            if len(candidates) == 1:
                e, d = candidates[0]
                prev = _unpack_note(e.note)
                assert prev is not None
                assigned[i] = e.record_id
                matched_entries.add(e.record_id)
                report.amended.append((e.record_id, NATURAL_KEY[d], prev[d], key[d]))
                e.natural_key_hash = key_hashes[i]
                e.content_hash = content_hashes[i]
                e.amendment_count += 1
                e.amended_from_key_hash = _hash(list(prev))
                e.last_seen_release = release
                e.last_modified_release = release
                e.last_modified_commit = commits.get(e.chunk_file, "")

    # Pass 3 -- allocate new ids. Each allocation is recorded against the line of the chunk it came
    # from, so the retirement pass below can answer "which id took over from the one I had".
    allocated_at: dict[tuple[str, int], list[tuple[str, tuple[str, ...]]]] = defaultdict(list)
    for i, a in enumerate(assigned):
        if a is not None:
            continue
        kh = key_hashes[i]
        first = seen_key.get(kh)
        if first is not None and assigned[first] is not None:
            assigned[i] = assigned[first]
            report.replicated.append((str(assigned[first]), str(rows[i].get("chunk.file") or "")))
            continue
        seen_key[kh] = i
        row = rows[i]
        chunk_file = str(row.get("chunk.file") or "")
        rid = registry._allocate()
        assigned[i] = rid
        commit = commits.get(chunk_file, "")
        registry._by_key[key_hashes[i]] = _Entry(
            record_id=rid, state=RecordState.ACTIVE,
            natural_key_hash=key_hashes[i], content_hash=content_hashes[i],
            chunk_file=chunk_file, chunk_row=int(row.get("chunk.row") or 0),
            first_seen_release=release, first_seen_commit=commit,
            last_seen_release=release, last_modified_release=release,
            last_modified_commit=commit,
            amendment_count=0, amended_from_key_hash="", replaced_by="",
            note=_pack_note(keys[i]),
        )
        report.added.append(rid)
        allocated_at[(chunk_file, int(row.get("chunk.row") or 0))].append((rid, keys[i]))

    # Re-key the registry (amendments changed natural_key_hash) and refresh the amendment notes.
    rekeyed: dict[str, _Entry] = {}
    for e in registry._by_key.values():
        rekeyed[e.natural_key_hash] = e
    registry._by_key = rekeyed
    # A different name from the `rid` allocated above: `assigned` holds `str | None` until every
    # slot is filled, and reusing the name widens it for the rest of the function.
    for i, settled in enumerate(assigned):
        e = registry.get(key_hashes[i])
        if e is not None and e.record_id == settled:
            e.note = _pack_note(keys[i])

    # Retire anything the registry still holds that this build did not see.
    # `added` is hoisted: rebuilding the set inside the loop makes this O(n^2), which on a 192k-row
    # corpus is the difference between a second and several minutes.
    added = set(report.added)
    retiring = [e for e in registry.active()
                if e.record_id not in matched_entries and e.record_id not in added]
    # How many are retiring from each line, so a line that lost two records links neither: with two
    # on each side there is no way to say which took over from which.
    retiring_at: Counter[tuple[str, int]] = Counter((e.chunk_file, e.chunk_row) for e in retiring)
    for e in retiring:
        e.state = RecordState.RETIRED
        e.replaced_by = _successor(e, allocated_at, retiring_at)
        report.retired.append(e.record_id)

    _backfill_successors(registry)

    out = records.with_columns(pl.Series("record_id", assigned, dtype=pl.String))
    return out, registry, report


def _backfill_successors(registry: IdentityRegistry) -> int:
    """Fill ``replaced_by`` on retirements that predate #693. Returns how many it filled.

    The evidence is still in the registry: an active record sitting at the line a retired one left,
    with the same receptor. The guards are :func:`_successor`'s, one step weaker in only one way -
    the successor is any active entry at that line rather than one allocated in this run, because
    the run that allocated it has already happened.

    Runs on every reconciliation and is idempotent: an entry that already has a pointer is skipped,
    and one that cannot be linked stays empty rather than being linked to something plausible.
    """
    active_at: dict[tuple[str, int], list[_Entry]] = defaultdict(list)
    orphan_at: Counter[tuple[str, int]] = Counter()
    for e in registry._by_key.values():
        at = (e.chunk_file, e.chunk_row)
        if e.state == RecordState.ACTIVE:
            active_at[at].append(e)
        elif not e.replaced_by:
            orphan_at[at] += 1

    filled = 0
    for e in registry._by_key.values():
        if e.state != RecordState.RETIRED or e.replaced_by:
            continue
        at = (e.chunk_file, e.chunk_row)
        successors = active_at.get(at, [])
        if len(successors) != 1 or orphan_at[at] != 1:
            continue
        previous, current = _unpack_note(e.note), _unpack_note(successors[0].note)
        if previous is None or current is None:
            continue
        if any(previous[i] != current[i] for i in RECEPTOR_FIELDS):
            continue
        e.replaced_by = successors[0].record_id
        filled += 1
    return filled


def _successor(entry: _Entry,
               allocated: dict[tuple[str, int], list[tuple[str, tuple[str, ...]]]],
               retiring: Counter[tuple[str, int]]) -> str:
    """The id that took over from a retiring one, or empty when nothing can be said. #693.

    The amendment pass keeps a `record_id` when the natural key moved in exactly one field, and
    refuses when two moved - correctly, because two fields moving is a different record by the rule
    the registry is built on. But the old id then had no pointer to the new one, which is the failure
    `ROADMAP.md` §10.4 exists to prevent: a reference that used to resolve returns nothing and a
    consumer cannot tell a correction from a deletion.

    The link is the same identity argument the `chunk.row` tie-break already makes - *the same line
    of the same file* - with two guards that make it a fact rather than a guess:

    * **one in, one out.** A line that retired two ids, or that had two allocated against it, links
      neither. There is no evidence for the pairing.
    * **the receptor is unchanged.** The two keys must agree on all six of
      :data:`RECEPTOR_FIELDS`. A new record at the line an old one left is a coincidence; the same
      TCR at that line, with only its annotation moved, is not.

    The case that motivated it is `#633`'s `RGPGRAFVTI` patch, where one line of
    `patches/antigen_epitope_species_gene.dict` moves `antigen.species` and `antigen.gene` together
    because both are one assertion about the peptide's source. The B16 `Plod1 -> Plod2` repair is
    49 more of the same shape: `antigen.gene` and `meta.subject.cohort` in one commit.
    """
    at = (entry.chunk_file, entry.chunk_row)
    candidates = allocated.get(at, [])
    if len(candidates) != 1 or retiring[at] != 1:
        return ""
    rid, key = candidates[0]
    previous = _unpack_note(entry.note)
    if previous is None:
        return ""
    if any(previous[i] != key[i] for i in RECEPTOR_FIELDS):
        return ""
    return rid


# The registry keeps the previous build's natural-key fields so an amendment can be *located*, not
# merely detected. Stored packed in `note` because the registry is a flat reviewable TSV and adding
# thirteen more columns would make its diffs unreadable.
# Must not be a character `str.splitlines()` or a CSV reader treats as a line break -- \x1c, \x1d and
# \x1e all are, and a note packed with one is silently truncated at its first field. \x1f is not.
_SEP = "\x1f"


def _pack_note(key: tuple[str, ...]) -> str:
    return _SEP.join(k.replace(_SEP, " ") for k in key)


def _unpack_note(note: str) -> tuple[str, ...] | None:
    if not note:
        return None
    parts = tuple(note.split(_SEP))
    if len(parts) == len(LEGACY_NATURAL_KEY):
        legacy = dict(zip(LEGACY_NATURAL_KEY, parts, strict=True))
        return tuple(legacy.get(c, "") for c in NATURAL_KEY)
    return parts if len(parts) == len(NATURAL_KEY) else None
