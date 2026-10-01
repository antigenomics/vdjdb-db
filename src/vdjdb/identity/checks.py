"""The identity invariants, as checks over a built directory.

``ROADMAP.md`` section 10.5 names seven. Six of them are properties of a build and live here; the
seventh, that a permuted chunk order changes no id, is a property of two builds and lives in
``tests/unit/test_identity_levels.py`` where it runs on a small fixture in seconds rather than
needing two full builds.

Each check returns rows rather than a boolean, so a failure says which ids and how many instead of
only that something is wrong.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from . import lifecycle as lc
from .levels import (
    CLONE,
    CLONOTYPE,
    EPITOPE,
    LEVELS,
    PMHC,
    Level,
    duplicate_ids,
    recomputed_mismatches,
)

#: The levels whose id is a hash of columns on the frame carrying it, so it can be recomputed in
#: place. ``clone`` is absent: its key is the pair of clonotype ids on *another* row, and
#: :func:`_clone_shape` checks it instead.
RECOMPUTABLE: tuple[tuple[str, Level], ...] = (
    ("chains", CLONOTYPE), ("records", PMHC), ("records", EPITOPE),
)


@dataclass(frozen=True, slots=True)
class Finding:
    """One failed invariant, with enough detail to act on."""

    invariant: int
    what: str
    rows: int
    sample: str = ""

    def __str__(self) -> str:
        tail = f"  e.g. {self.sample}" if self.sample else ""
        return f"invariant {self.invariant}: {self.what}, {self.rows:,} rows{tail}"


def _sample(frame: pl.DataFrame, n: int = 3) -> str:
    if frame.is_empty():
        return ""
    return "; ".join(str(r) for r in frame.head(n).rows())


def _distinct(tables: dict[str, pl.DataFrame]) -> list[Finding]:
    """Invariant 1: no id is shared by two distinct keys."""
    out = []
    for table, level in RECOMPUTABLE:
        frame = tables.get(table)
        if frame is None or level.column not in frame.columns:
            continue
        # `species` is on `records`, so the clonotype key is only complete after the join below.
        key = [c for c in level.key if c in frame.columns]
        dup = duplicate_ids(frame.filter(pl.col(level.column) != ""), level.column, key)
        if not dup.is_empty():
            out.append(Finding(1, f"{level.name}_id collides on {len(key)} key columns",
                               dup.height, _sample(dup)))
    return out


def _recomputed(tables: dict[str, pl.DataFrame]) -> list[Finding]:
    """Invariant 2: every stored id equals the hash of the key on its own row."""
    out = []
    records = tables.get("records")
    for table, level in RECOMPUTABLE:
        frame = tables.get(table)
        if frame is None or level.column not in frame.columns:
            continue
        if table == "chains" and records is not None and "species" not in frame.columns:
            frame = frame.join(records.select("record_id", "species"), on="record_id",
                               how="left", maintain_order="left")
        if any(c not in frame.columns for c in level.key):
            continue
        bad = recomputed_mismatches(frame, level)
        if not bad.is_empty():
            out.append(Finding(2, f"{level.name}_id is not the hash of its own key",
                               bad.height, _sample(bad)))
    return out


def _clone_shape(tables: dict[str, pl.DataFrame]) -> list[Finding]:
    """Invariant 6, the clone half: a clone covers exactly two chains of different ``gene``.

    Checked per record, not per id: one clone is reported by many records, and the property that
    matters is that each record's pair is a pair.
    """
    chains = tables.get("chains")
    if chains is None or CLONE.column not in chains.columns:
        return []
    # `fill_null` first: a null would otherwise make the comparison null and drop the row silently,
    # and a null in a string column is itself a rule-6 violation that a different test catches.
    shape = (chains.filter(pl.col(CLONE.column).fill_null("") != "")
                   .group_by("record_id", CLONE.column)
                   .agg(pl.len().alias("n"), pl.col("gene").n_unique().alias("genes")))
    bad = shape.filter((pl.col("n") != 2) | (pl.col("genes") != 2))
    if bad.is_empty():
        return []
    return [Finding(6, "a clone_id does not cover exactly two chains of different gene",
                    bad.height, _sample(bad))]


def _referential(tables: dict[str, pl.DataFrame]) -> list[Finding]:
    """Invariant 6, the rest: every id referenced somewhere resolves where it should."""
    out = []
    records, chains, evidence = (tables.get(n) for n in ("records", "chains", "evidence"))
    if records is not None and chains is not None:
        orphan = chains.join(records.select("record_id"), on="record_id", how="anti")
        if not orphan.is_empty():
            out.append(Finding(6, "a chain names a record_id that records does not have",
                               orphan.height, _sample(orphan.select("record_id").unique())))
    if records is not None:
        for level in (PMHC, EPITOPE):
            if level.column not in records.columns:
                continue
            blank = records.filter((pl.col(level.column) == "")
                                   | pl.col(level.column).is_null())
            if not blank.is_empty():
                out.append(Finding(6, f"{level.column} is empty on a record",
                                   blank.height, _sample(blank.select("record_id"))))
    if evidence is not None and chains is not None and "gene" in evidence.columns:
        orphan = evidence.join(chains.select("record_id", "gene"), on=["record_id", "gene"],
                               how="anti")
        if not orphan.is_empty():
            out.append(Finding(6, "an evidence row names a chain that does not exist",
                               orphan.height, _sample(orphan.select("record_id", "gene").unique())))
    return out


def _against_previous(tables: dict[str, pl.DataFrame], previous: pl.DataFrame) -> list[Finding]:
    """Invariant 5: an id that comes back is the same id, and no prefix is unknown.

    A derived id returning is legitimate, because the id *is* the hash of its key, so the same key
    returning has to produce the same id. What would be reuse is an id resolving to a different key,
    and invariant 2 is what rules that out. So this check is the cheap complement: every id in the
    build whose prefix belongs to no level, which is how a hand-edited fixture drifts.
    """
    now = lc.present(tables)
    out = []
    unknown = lc.unknown_prefixes(now)
    if unknown:
        out.append(Finding(5, "an id carries a prefix belonging to no level",
                           len(unknown), "; ".join(unknown[:3])))
    # A level that was published and is now entirely absent is not an error, but it is never
    # accidental either, so it is reported as a finding for a human to accept.
    for level in LEVELS:
        was = previous.filter((pl.col("level") == level.name)
                              & (pl.col("state") == lc.ACTIVE)).height
        is_ = now.filter(pl.col("level") == level.name).height
        if was and not is_:
            out.append(Finding(5, f"level {level.name} had {was:,} active ids and now has none", was))
    return out


def _tcr_hash_unchanged(tables: dict[str, pl.DataFrame],
                        previous_chains: pl.DataFrame) -> list[Finding]:
    """Invariant 7: ``TCR_hash`` is byte-identical wherever the chain key is unchanged.

    The legacy structure link is stored as curated and never recomputed, so a record whose key did
    not move must carry the same hash it shipped with. Keyed on ``clonotype_id`` rather than on
    ``record_id``, because the structure is a property of the receptor.
    """
    chains = tables.get("chains")
    if chains is None or "TCR_hash" not in chains.columns:
        return []
    if "TCR_hash" not in previous_chains.columns or CLONOTYPE.column not in previous_chains.columns:
        return []
    before = (previous_chains.select(CLONOTYPE.column, "TCR_hash")
                             .filter(pl.col("TCR_hash") != "").unique())
    after = (chains.select(CLONOTYPE.column, "TCR_hash")
                   .filter(pl.col("TCR_hash") != "").unique())
    moved = (before.join(after, on=CLONOTYPE.column, how="inner", suffix="_now")
                   .filter(pl.col("TCR_hash") != pl.col("TCR_hash_now")))
    if moved.is_empty():
        return []
    return [Finding(7, "TCR_hash moved on a clonotype whose key did not", moved.height,
                    _sample(moved))]


def check(tables: dict[str, pl.DataFrame], *,
          previous_lifecycle: pl.DataFrame | None = None,
          previous_chains: pl.DataFrame | None = None) -> list[Finding]:
    """Every invariant this directory can answer, worst first by row count.

    The checks needing history are skipped rather than failed when it is absent: a fork and a first
    build have none, and they are expected to pass.
    """
    found = (_distinct(tables) + _recomputed(tables) + _clone_shape(tables)
             + _referential(tables))
    if previous_lifecycle is not None and not previous_lifecycle.is_empty():
        found += _against_previous(tables, previous_lifecycle)
    if previous_chains is not None:
        found += _tcr_hash_unchanged(tables, previous_chains)
    return sorted(found, key=lambda f: (-f.rows, f.invariant))


def read_tables(directory: Path) -> dict[str, pl.DataFrame]:
    """Whatever of the definitive tables the directory holds. Parquet, which is the shipped form.

    Unlike :func:`vdjdb.emit.vdjdb3.read_tables` this tolerates a partial directory, because the
    identity checks are meant to run on whatever a build got as far as writing. The per-table read
    is that module's, so the shipped ``underscore_case`` header is translated in one place.
    """
    from ..emit.vdjdb3 import read_table

    out = {}
    for name in ("records", "chains", "evidence", "epitopes", "restriction"):
        if (Path(directory) / f"{name}.parquet").exists():
            out[name] = read_table(Path(directory), name)
    return out
