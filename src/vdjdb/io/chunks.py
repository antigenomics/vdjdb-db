"""Read ``chunks/`` into one polars frame.

``chunks/`` **is** the data: one file per publication, 230 files, ~203k rows. Everything else in the
repository is machinery for validating, assembling and publishing it.

Three properties this reader guarantees that the pandas one did not:

* **Determinism.** ``runBuidDatabase.py`` iterates ``os.listdir``, i.e. readdir order, which is
  filesystem- and host-dependent. That is the single reason the 2026-06-03 release cannot be
  reproduced byte-for-byte by anyone, including the pipeline that produced it. This reader sorts.
* **One missing marker.** Every cell is a string and every absent cell is ``""``. The pandas path's
  ``None`` / ``NaN`` / ``""`` three-way ambiguity is the source of more than one shipped bug.
* **Provenance.** Every row carries ``chunk.file`` and ``chunk.row``, so a record can be traced back
  to the line a curator wrote. The identity registry keys amendments on it.

Deduplication is **within a chunk**, and that is not a performance choice. A chunk is one paper, so
two matching rows in two chunks are two papers reporting the same receptor independently -- the
strongest evidence the database carries, and the signal motif clustering is tuned against
(ROADMAP section 11.1). Collapsing them would delete it. Measured: 19 such pairs.
"""
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ..config import Paths
from ..schema import ALL_COLUMNS, CHUNK_DEDUP_KEY, KEPT_CURATION_COLUMNS

#: Provenance columns this reader adds. Not part of any chunk.
PROVENANCE: tuple[str, ...] = ("chunk.file", "chunk.row")

#: Columns read from a chunk when present: the 31 required plus the five kept for debugging.
READABLE: tuple[str, ...] = ALL_COLUMNS + KEPT_CURATION_COLUMNS


def chunk_files(directory: Path | None = None) -> list[Path]:
    """Every chunk file, **sorted**. Hidden files are skipped, as the legacy build skipped them."""
    d = directory or Paths.discover().chunks
    return sorted(p for p in d.iterdir()
                  if p.suffix in {".txt", ".tsv"} and not p.name.startswith("."))


def _normalise_header(name: str) -> str:
    """Strip the BOM and surrounding whitespace, and lower-case the one capitalised column.

    ``vandesandt-etal-2019-11-04.txt`` writes ``Comment``. Nineteen distinct header rows exist
    across 230 files; the rest of the normalisation is a one-shot migration (#497, phase 3), not a
    read-time transform, because silently accepting a malformed header is how they accumulated.
    """
    n = name.lstrip("﻿").strip()
    return "comment" if n == "Comment" else n


def read_chunk(path: Path) -> pl.DataFrame:
    """One chunk as an all-string frame with provenance, missing optional columns filled with ``""``.

    ``chunk.row`` is 1-based and counts data rows, so it is the line number a curator sees in a
    spreadsheet minus the header.
    """
    df = pl.read_csv(
        path,
        separator="\t",
        quote_char=None,          # no field is ever quoted; see tests/release
        has_header=True,
        infer_schema_length=0,    # every cell a string: "0" and "0.0" are different cells
        encoding="utf8-lossy",    # matches the legacy encoding_errors="ignore"
        truncate_ragged_lines=True,
        null_values=[],
    )
    # A duplicate column name would make the selection below ambiguous. polars silently renames
    # the second one (`species_duplicated_0`), so the raw header is the only place to catch it.
    with path.open("rb") as fh:
        raw_header = [_normalise_header(c)
                      for c in fh.readline().decode("utf-8", "replace").rstrip("\r\n").split("\t")]
    if len(set(raw_header)) != len(raw_header):
        dupes = sorted({c for c in raw_header if raw_header.count(c) > 1})
        raise ValueError(f"{path.name}: duplicate columns {dupes}")

    df = df.rename({c: _normalise_header(c) for c in df.columns})

    missing = [c for c in ALL_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing required columns {missing}")

    present = [c for c in READABLE if c in df.columns]
    return (
        df.select(present)
        # Strip surrounding whitespace, not only the CR a CRLF file leaves behind. It is never
        # meaningful in a TSV cell and it silently forks a value in two: measured, 758 record-cells
        # across 8 columns, including `tetramer-sort ` appearing beside `tetramer-sort` (103 records)
        # and `Nucleocapsid ` (171). `strip_chars()` with no argument also takes the non-breaking
        # space that one J-gene call carried.
        .with_columns(pl.col(present).cast(pl.Utf8).fill_null("").str.strip_chars())
        .with_columns(
            *(pl.lit("").alias(c) for c in READABLE if c not in present),
            pl.lit(path.name).alias("chunk.file"),
            (pl.int_range(1, pl.len() + 1, dtype=pl.UInt32)).alias("chunk.row"),
        )
        .select(*READABLE, *PROVENANCE)
    )


def dedup(df: pl.DataFrame) -> pl.DataFrame:
    """Per-chunk deduplication on :data:`CHUNK_DEDUP_KEY`, keeping the first occurrence.

    Per-chunk, not global: the released ``vdjdb_full.txt`` is the per-chunk result, and switching to
    global would silently drop 19 rows.
    """
    return df.unique(subset=["chunk.file", *CHUNK_DEDUP_KEY], keep="first", maintain_order=True)


def read_chunks(paths: Iterable[Path] | None = None, *, deduplicate: bool = True) -> pl.DataFrame:
    """Every chunk, concatenated in sorted filename order.

    Measured: 230 files, 203,308 raw rows, 192,753 after per-chunk deduplication -- exactly the
    released row count -- in 0.4 s.
    """
    files = list(paths) if paths is not None else chunk_files()
    if not files:
        raise ValueError("no chunk files found")
    df = pl.concat([read_chunk(p) for p in files], how="vertical")
    return dedup(df) if deduplicate else df


def independently_reported(df: pl.DataFrame) -> pl.DataFrame:
    """Records reported by more than one chunk -- that is, by more than one paper.

    **Not duplicates.** A chunk is one publication, so the same receptor against the same epitope
    appearing in two chunks is independent replication. This is evidence, and phase 11 fits the
    motif clustering against it.
    """
    counts = df.group_by(CHUNK_DEDUP_KEY).agg(
        pl.col("chunk.file").n_unique().alias("chunks"),
        pl.col("chunk.file").unique().sort().str.join(",").alias("chunk.files"),
    )
    return counts.filter(pl.col("chunks") > 1).sort("chunk.files")
