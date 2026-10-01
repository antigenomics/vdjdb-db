"""Read ``chunks/`` into one polars frame.

``chunks/`` is the data: one file per publication, 230 files, ~203k rows. Everything else in the
repository is machinery for validating, assembling and publishing it.

Three properties this reader guarantees that the pandas one did not:

* **Determinism.** ``runBuidDatabase.py`` iterates ``os.listdir``, i.e. readdir order, which is
  filesystem- and host-dependent. That is the single reason the 2026-06-03 release cannot be
  reproduced byte-for-byte by anyone, including the pipeline that produced it. This reader sorts.
* **One missing marker.** Every cell is a string and every absent cell is ``""``. The pandas path's
  ``None`` / ``NaN`` / ``""`` three-way ambiguity is the source of more than one shipped bug.
* **Provenance.** Every row has ``chunk.file`` and ``chunk.row``, so a record can be traced back
  to the line a curator wrote. The identity registry keys amendments on it.

Deduplication is within a chunk, not globally. A chunk is one paper, so two matching rows in two
chunks are two papers reporting the same receptor independently -- the strongest evidence the
database holds, and the signal motif clustering is tuned against (ROADMAP section 11.1). Collapsing
them would delete it. Measured: 19 such pairs.
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
    """Every chunk file, sorted. Hidden files are skipped, as the legacy build skipped them."""
    d = directory or Paths.discover().chunks
    return sorted(p for p in d.iterdir()
                  if p.suffix in {".txt", ".tsv"} and not p.name.startswith("."))


def _normalise_header(name: str) -> str:
    """Strip the BOM and surrounding whitespace, and lower-case the one capitalised column.

    ``vandesandt-etal-2019-11-04.tsv`` writes ``Comment``. Nineteen distinct header rows exist
    across 230 files; the rest of the normalisation is a one-shot migration (#497, phase 3), not a
    read-time transform, because accepting a malformed header without an error is how they
    accumulated.
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
    # A duplicate column name would make the selection below ambiguous. polars renames the second
    # one (`species_duplicated_0`) without an error, so the raw header is where to catch it.
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
        # meaningful in a TSV cell and it forks a value in two: measured, 758 record-cells across
        # 8 columns, including `tetramer-sort ` appearing beside `tetramer-sort` (103 records) and
        # `Nucleocapsid ` (171). `strip_chars()` with no argument also takes the non-breaking space
        # that one J-gene call had.
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
    global would drop 19 rows.
    """
    return df.unique(subset=["chunk.file", *CHUNK_DEDUP_KEY], keep="first", maintain_order=True)


#: Columns that describe the **curation** rather than the record, so a difference in one is not a
#: difference between two reports. `chunk.id` is a per-chunk row serial, so two chunks curating one
#: clone always disagree on it; before it was excluded, all 19 of #390's groups read as conflicting
#: for that reason alone.
CURATION_ONLY: frozenset[str] = frozenset({"chunk.file", "chunk.row", "chunk.id",
                                           "submitter", "comment"})

#: Chunks that are not one publication's own report. A structure chunk and an aggregate both carry
#: rows whose paper has its own chunk, so where a merge has to pick a base row, it prefers the
#: paper's. Matched by name; anything else is a submitted chunk.
DERIVED_CHUNKS: tuple[str, ...] = ("PDB_Database.tsv", "small_datasets_")


def _base_first(files: list[str]) -> list[str]:
    """``files`` with the submitted paper chunk first - trust the submission first.

    Deterministic where neither is derived, which the two `goncharov-*` pairs are: sorted by name
    after the derived ones are pushed back, never by row order (hard rule 7).
    """
    return sorted(files, key=lambda f: (f.startswith(DERIVED_CHUNKS), f))


def merge_repeated_references(df: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Collapse rows that are one publication reporting one clone twice, in two chunk files (#390).

    ``CHUNK_DEDUP_KEY`` contains ``reference.id``, so a group of it spanning two chunk files is one
    paper curated twice - and two rows of one paper are not two independent reports, which is the
    whole basis for deduplicating within a chunk in the first place. **The chunk was a proxy for the
    publication**; where the two come apart, the publication is what counts.

    Measured on the corpus: **19 groups over 38 rows**, every one a pair. 18 agree on every column
    that describes the record and are merged, filling the base row's blanks from the other -
    ``method.verification`` in 10 of the 18 and ``meta.epitope.id`` in 3. **1 is left alone**:
    `PDB_Database.tsv` and `PMID_34433824.tsv` give one clone `structural` and `tetramer-sort`, which
    is one paper reporting a solved complex *and* the sort that found it - two observations, not a
    duplicate.

    A group is merged only when every column both rows carry non-blank agrees, ``CURATION_ONLY``
    excluded. Nothing is overwritten: the base row keeps every value it has.
    """
    report_schema = {"reference.id": pl.String, "chunk.files": pl.String, "verdict": pl.String,
                     "detail": pl.String, "rows": pl.Int64}
    spanning = (df.group_by(CHUNK_DEDUP_KEY)
                  .agg(pl.col("chunk.file").n_unique().alias("__files"))
                  .filter(pl.col("__files") > 1)
                  .drop("__files"))
    if spanning.is_empty():
        return df, pl.DataFrame(schema=report_schema)

    compared = [c for c in df.columns if c not in CHUNK_DEDUP_KEY and c not in CURATION_ONLY]
    groups = df.join(spanning, on=CHUNK_DEDUP_KEY, how="inner")
    rows, drop, patch = [], [], []
    for key, sub in groups.group_by(CHUNK_DEDUP_KEY, maintain_order=True):
        files = _base_first(sub["chunk.file"].to_list())
        disagree, fill = {}, {}
        for column in compared:
            values = {v for v in sub[column].to_list() if v not in (None, "")}
            if len(values) > 1:
                disagree[column] = " | ".join(sorted(values))
            elif values:
                fill[column] = next(iter(values))
        reference = key[CHUNK_DEDUP_KEY.index("reference.id")]
        if disagree:
            rows.append({"reference.id": reference, "chunk.files": ",".join(files),
                         "verdict": "kept", "rows": sub.height,
                         "detail": "; ".join(f"{c}: {v}" for c, v in sorted(disagree.items()))})
            continue
        base = sub.filter(pl.col("chunk.file") == files[0]).row(0, named=True)
        rows.append({"reference.id": reference, "chunk.files": ",".join(files),
                     "verdict": "merged", "rows": sub.height,
                     "detail": ("filled " + ", ".join(sorted(c for c in fill if not base[c]))
                                if any(not base[c] for c in fill) else "identical")})
        drop += [(f, r) for f, r in zip(sub["chunk.file"], sub["chunk.row"], strict=True)
                 if not (f == base["chunk.file"] and r == base["chunk.row"])]
        patch.append({**base, **{c: v for c, v in fill.items() if not base[c]}})

    if drop:
        dropped = pl.DataFrame(drop, schema={"chunk.file": pl.String, "chunk.row": df["chunk.row"].dtype},
                               orient="row")
        kept = df.join(dropped, on=["chunk.file", "chunk.row"], how="anti")
        patched = pl.DataFrame(patch, schema=df.schema)
        df = (pl.concat([kept.join(patched.select("chunk.file", "chunk.row"),
                                   on=["chunk.file", "chunk.row"], how="anti"), patched])
                .sort("chunk.file", "chunk.row"))
    report = pl.DataFrame(rows, schema=report_schema)
    return df, report.sort("verdict", "reference.id", "chunk.files")


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

    Not duplicates. A chunk is one publication, so the same receptor against the same epitope
    appearing in two chunks is independent replication. This is evidence, and phase 11 fits the
    motif clustering against it.
    """
    counts = df.group_by(CHUNK_DEDUP_KEY).agg(
        pl.col("chunk.file").n_unique().alias("chunks"),
        pl.col("chunk.file").unique().sort().str.join(",").alias("chunk.files"),
    )
    return counts.filter(pl.col("chunks") > 1).sort("chunk.files")
