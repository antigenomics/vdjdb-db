"""Stable VDJdb record identifiers.

The problem this solves: a VDJdb record has no identifier today. It is a row in a chunk file, and
the only thing that names it is its own content. So when a curator fixes a typo in a CDR3, the old
record silently disappears and a new one silently appears, and nothing downstream can tell that from
a genuine deletion plus a genuine addition. Citations to vdjdb.com, structure links keyed on
``TCR_hash``, and evidence accumulated across releases all break invisibly.

Identity here is therefore **two-level**:

``record_id``
    Opaque, assigned once, never reused, and *stable across content changes*. This is what external
    references point at. Format ``VDJDB<10 digits>``, allocated monotonically.

``content_hash``
    sha256 over the record's canonical content. Changes whenever anything about the record changes.
    This is what *detects* a change; it is never an identifier.

Plus a lifecycle (``active`` / ``amended`` / ``retired``) and provenance back to the chunk file, the
row within it, and the commit that last touched it.

**Matching cascade.** On each build, records from ``chunks/`` are reconciled against the committed
registry in three passes, most specific first:

1. **Exact** natural key -> the same ``record_id``, unchanged.
2. **Amendment**: within the same chunk and the same ``reference.id``, a record whose natural key
   differs from an unmatched registry entry in **exactly one field** is treated as that entry,
   amended. This is what makes a typo fix traceable instead of a delete-plus-insert. Ambiguous
   matches (two equally good candidates) are refused and fall through -- a wrong link is worse than
   a new id.
3. **New**: anything still unmatched gets a fresh id. Registry entries still unmatched afterwards
   are ``retired``, carrying the release in which they were last seen.

The registry is a committed TSV sorted by ``record_id`` so that a curation PR shows exactly which
records were added, amended or retired, as a reviewable diff.
"""
from __future__ import annotations

import hashlib
from collections import defaultdict
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

import polars as pl

from ..schema import ALL_COLUMNS, CHUNK_DEDUP_KEY

#: The fields that identify a record. A change in any of these is an *amendment*; a change anywhere
#: else is re-annotation that keeps the same record.
#:
#: ``CHUNK_DEDUP_KEY`` **plus the chunk**. A chunk is one paper, so two rows in two chunks are
#: independent reports rather than one record seen twice, however identical their fields (CLAUDE.md,
#: the data model). Deduplication is within a chunk; identity is per curated line.
#:
#: Two earlier versions were wrong in opposite directions. A narrower key that stopped at
#: ``reference.id`` collided on **20,769 of 192,753** records -- one paper reporting the same TCR
#: against the same epitope in several donors is several records. Dropping ``chunk.file`` merged
#: **19** pairs that are two papers' independent reports, which is exactly the signal phase 11
#: tunes against.
NATURAL_KEY: tuple[str, ...] = (*CHUNK_DEDUP_KEY, "chunk.file")

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
    "note",
)


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

    Changes whenever anything about the record changes -- including method and meta, which the
    natural key ignores. That is the point: the natural key answers "is this the same record", the
    content hash answers "has it changed".
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
        # Empty string is the only missing marker (CLAUDE.md), and the polars kwarg that spells
        # that has been renamed across versions -- normalise after the read instead.
        df = pl.read_csv(path, separator="\t", infer_schema=False, quote_char=None).fill_null("")
        missing = set(REGISTRY_COLUMNS) - set(df.columns)
        if missing:
            raise ValueError(f"{path} is missing registry columns: {sorted(missing)}")
        entries = [
            _Entry(
                record_id=r["record_id"], state=r["state"],
                natural_key_hash=r["natural_key_hash"], content_hash=r["content_hash"],
                chunk_file=r["chunk_file"], chunk_row=int(r["chunk_row"] or 0),
                first_seen_release=r["first_seen_release"], first_seen_commit=r["first_seen_commit"],
                last_seen_release=r["last_seen_release"],
                last_modified_release=r["last_modified_release"],
                last_modified_commit=r["last_modified_commit"],
                amendment_count=int(r["amendment_count"] or 0),
                amended_from_key_hash=r["amended_from_key_hash"], note=r["note"],
            )
            for r in df.iter_rows(named=True)
        ]
        return cls(entries)

    def save(self, path: Path) -> None:
        rows = sorted(self._by_key.values(), key=lambda e: parse_id(e.record_id))
        df = pl.DataFrame(
            {c: [str(getattr(e, c)) for e in rows] for c in REGISTRY_COLUMNS}
        ) if rows else pl.DataFrame({c: pl.Series(c, [], pl.String) for c in REGISTRY_COLUMNS})
        path.parent.mkdir(parents=True, exist_ok=True)
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

    ``records`` must carry ``ALL_COLUMNS`` plus ``chunk.file`` and ``chunk.row``.
    ``commits`` maps a chunk filename to the commit that last touched it.
    """
    commits = commits or {}
    report = ReconcileReport()

    rows = records.to_dicts()
    keys = [_key_fields(r) for r in rows]
    key_hashes = [natural_key(r) for r in rows]
    content_hashes = [canonical_content_hash(r) for r in rows]

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
        ref_idx = NATURAL_KEY.index("reference.id")
        leftovers = [e for e in registry.active() if e.record_id not in matched_entries]
        # Bucket candidates by (chunk_file, reference.id) so the scan stays local.
        by_bucket: dict[tuple[str, str], list[_Entry]] = defaultdict(list)
        for e in leftovers:
            prev = _unpack_note(e.note)
            if prev is not None:
                by_bucket[(e.chunk_file, prev[ref_idx])].append(e)

        # The registry does not store the raw key, only its hash, so amendment matching needs the
        # previous build's key fields. They are carried in `note` as the reference id plus the
        # packed key; see `_pack_note`. Entries written by older versions simply do not match.
        for i in unmatched_rows:
            row, key = rows[i], keys[i]
            bucket = by_bucket.get((str(row.get("chunk.file") or ""), key[ref_idx]), [])
            candidates: list[tuple[_Entry, int]] = []
            for e in bucket:
                if e.record_id in matched_entries:
                    continue
                prev = _unpack_note(e.note)
                if prev is None:
                    continue
                d = _single_field_difference(key, prev)
                if d is not None:
                    candidates.append((e, d))
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

    # Pass 3 -- allocate new ids.
    for i, a in enumerate(assigned):
        if a is not None:
            continue
        kh = key_hashes[i]
        first = seen_key.get(kh)
        if first is not None and assigned[first] is not None:
            assigned[i] = assigned[first]
            report.duplicated.append((str(assigned[first]), str(rows[i].get("chunk.file") or "")))
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
            amendment_count=0, amended_from_key_hash="", note=_pack_note(keys[i]),
        )
        report.added.append(rid)

    # Re-key the registry (amendments changed natural_key_hash) and refresh the amendment notes.
    rekeyed: dict[str, _Entry] = {}
    for e in registry._by_key.values():
        rekeyed[e.natural_key_hash] = e
    registry._by_key = rekeyed
    for i, rid in enumerate(assigned):
        e = registry.get(key_hashes[i])
        if e is not None and e.record_id == rid:
            e.note = _pack_note(keys[i])

    # Retire anything the registry still holds that this build did not see.
    # `added` is hoisted: rebuilding the set inside the loop makes this O(n^2), which on a 192k-row
    # corpus is the difference between a second and several minutes.
    added = set(report.added)
    for e in registry.active():
        if e.record_id not in matched_entries and e.record_id not in added:
            e.state = RecordState.RETIRED
            report.retired.append(e.record_id)

    out = records.with_columns(pl.Series("record_id", assigned, dtype=pl.String))
    return out, registry, report


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
    return parts if len(parts) == len(NATURAL_KEY) else None
