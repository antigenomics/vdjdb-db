"""Size-limited generated text tables and readers for either on-disk form."""
from __future__ import annotations

import gzip
import os
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO

import polars as pl

MAX_TEXT_BYTES = 1_000_000


def resolve_path(path: Path | str) -> Path:
    """Accept a logical plain filename when only its gzip form exists."""
    path = Path(path)
    compressed = path.with_name(path.name + ".gz")
    return compressed if not path.exists() and compressed.exists() else path


def open_text(path: Path | str) -> TextIO:
    path = resolve_path(path)
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_csv(path, **kwargs) -> pl.DataFrame:
    """Polars reads gzip natively; preserve support for bytes and open streams."""
    if isinstance(path, (Path, str)):
        path = resolve_path(path)
    return pl.read_csv(path, **kwargs)


@contextmanager
def table_writer(path: Path | str) -> Iterator[TextIO]:
    """Publish atomically, compressing TSV/TXT above the limit and explicit gzip targets.

    Temporary files have no table suffix. A failed write leaves the previous output intact.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: list[Path] = []
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=path.parent, suffix=".tmp", delete=False) as stream:
            plain = Path(stream.name)
            temporary.append(plain)
            yield stream
        explicit = path.suffix == ".gz"
        compress = explicit or (path.suffix in {".tsv", ".txt"}
                                and plain.stat().st_size > MAX_TEXT_BYTES)
        target = path if explicit or not compress else path.with_name(path.name + ".gz")
        ready = plain
        if compress:
            with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".tmp", delete=False) as raw:
                ready = Path(raw.name)
                temporary.append(ready)
                with plain.open("rb") as source, gzip.GzipFile(
                        filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
                    shutil.copyfileobj(source, zipped)
        os.replace(ready, target)
        alternate = (path.with_suffix("") if explicit else
                     path if compress else path.with_name(path.name + ".gz"))
        alternate.unlink(missing_ok=True)
    finally:
        for tmp in temporary:
            tmp.unlink(missing_ok=True)


def write_csv(frame: pl.DataFrame, path: Path | str, **kwargs) -> Path:
    with table_writer(path) as stream:
        frame.write_csv(stream.buffer, **kwargs)
    return resolve_path(path)


def write_text(path: Path | str, content: str) -> Path:
    with table_writer(path) as stream:
        stream.write(content)
    return resolve_path(path)
