"""The new VDJdb format: the definitive tables as they ship, plus the joined convenience view.

Three fact tables and one view, parquet with a TSV projection of each (`docs/outputs.md` §3):

=====================  =========================================================================
``records``            one row per curated record, PK ``record_id``
``chains``             one row per TCR chain, PK ``(record_id, gene)``
``evidence``           one row per piece of support, PK ``(record_id, evidence_id)``
``vdjdb``              ``records`` |><| ``chains`` |><| ``evidence``, the denormalised view
=====================  =========================================================================

The first three are written, not built: :mod:`vdjdb.assemble.tables` produced them and this module
only chooses a file format. ``vdjdb`` is the one thing derived here, and it is derived on every
build, so an edit to it is not an edit to the database.

Unlike :mod:`vdjdb.emit.legacy`, nothing here reproduces a byte layout, so the writers are plain:
default quoting rather than ``quote_style="never"`` (a field containing a tab must survive, not
corrupt the row), and no JSON blobs at all.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import polars as pl

from ..schema import (
    CHAIN_COLUMNS,
    EPITOPE_COLUMNS,
    EVIDENCE_TABLE_COLUMNS,
    RECORD_COLUMNS,
    RESTRICTION_COLUMNS,
    TIDY_NAMES,
    schema_json,
)

#: ``evidence_type`` -> the boolean column it becomes in the joined view. Production ``vdjdb-web``
#: has served five ``evidence.*`` columns for years that nothing in this repo produced; this is
#: where they start being produced, so the names are theirs.
EVIDENCE_VIEW: dict[str, str] = {
    "independent_study": "evidence.validation.independent",
    "motif_tcrnet": "evidence.motif.tcrnet",
    "motif_tcremp": "evidence.motif.tcremp",
    "structure_native": "evidence.structure.native",
    "structure_model": "evidence.structure.model",
}

#: Declared in full, so the view's schema does not change shape as producers land. Everything
#: without a producer yet is ``false`` -- "no evidence of this kind", not a missing column.
#: ``same.study`` has no producer at all: method-level self-validation is phase 9's.
VIEW_EVIDENCE_COLUMNS: tuple[str, ...] = (*EVIDENCE_VIEW.values(), "evidence.validation.same.study")

#: :data:`vdjdb.schema.TIDY_NAMES` extended with the six columns above. They are `vdjdb-web`'s own
#: names and the registry deliberately does not declare them, so the extension belongs here beside
#: the declaration rather than in the registry. Anything outside both still raises on write, which
#: is what catches a genuinely undeclared column.
_TIDY: dict[str, str] = {**TIDY_NAMES,
                         **{c: c.replace(".", "_") for c in VIEW_EVIDENCE_COLUMNS}}
_FROM_TIDY: dict[str, str] = {v: k for k, v in _TIDY.items()}

_TABLE_ORDER: dict[str, tuple[str, ...]] = {
    "records": RECORD_COLUMNS, "chains": CHAIN_COLUMNS, "evidence": EVIDENCE_TABLE_COLUMNS,
    "epitopes": EPITOPE_COLUMNS, "restriction": RESTRICTION_COLUMNS,
}


def joined(tables: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """``records`` joined to ``chains``, with each evidence type pivoted to a boolean column.

    One row per chain -- the level ``vdjdb.txt`` is written at -- so this is the table a user who
    wants one flat thing should read, and it is the one table here that is purely derived.
    """
    records, chains, evidence = tables["records"], tables["chains"], tables["evidence"]
    out = chains.join(records, on="record_id", how="left")

    ev = evidence.with_columns(
        pl.col("evidence_type").replace_strict(EVIDENCE_VIEW, default=None).alias("__col")
    ).drop_nulls("__col")
    if ev.select(pl.col("gene").fill_null("").eq("").any()).item():
        # Record-level evidence (a structure covers the receptor, not one chain) needs a second
        # join on `record_id` alone. Nothing produces it yet; raise rather than drop the rows.
        raise NotImplementedError(
            "record-level evidence rows (empty `gene`) need a join on record_id; "
            "add it here when the first producer lands")
    if not ev.is_empty():
        wide = (ev.select("record_id", "gene", "__col").unique()
                .with_columns(pl.lit(True).alias("__v"))
                .pivot(on="__col", index=["record_id", "gene"], values="__v",
                       aggregate_function="first"))
        out = out.join(wide, on=["record_id", "gene"], how="left")

    return out.with_columns(
        *(pl.col(c).fill_null(False) if c in out.columns else pl.lit(False).alias(c)
          for c in VIEW_EVIDENCE_COLUMNS)
    ).select(*CHAIN_COLUMNS, *(c for c in RECORD_COLUMNS if c != "record_id"),
             *VIEW_EVIDENCE_COLUMNS).sort("record_id", "gene")


def read_table(d: Path, name: str, columns: Sequence[str] | None = None) -> pl.DataFrame:
    """One tidy table from a build directory, under the **internal** dotted column names.

    The files ship ``underscore_case`` (:data:`vdjdb.schema.TIDY_NAMES`) and everything inside
    ``assemble/`` speaks the dotted names, so the rename lives here and at :func:`write_all` and
    nowhere else. Every reader of a build directory goes through this function for that reason: a
    second ``pl.read_parquet`` on these files is a second place the two conventions meet.

    ``columns`` is given in internal names and translated, so a caller still asks for what it means.
    """
    want = [_TIDY[c] for c in columns] if columns is not None else None
    return to_internal(pl.read_parquet(d / f"{name}.parquet", columns=want))


def to_internal(frame: pl.DataFrame) -> pl.DataFrame:
    """A tidy-table frame read off disk, put back under the internal dotted column names.

    :func:`read_table` is the way to read one of these tables. This is for the caller that has to
    read the TSV projection instead -- the legacy parity check does, because it is asserting on what
    shipped as text -- so that it does not need a second copy of the mapping.
    """
    return frame.rename({c: _FROM_TIDY[c] for c in frame.columns if c in _FROM_TIDY})


def read_tables(d: Path) -> dict[str, pl.DataFrame]:
    """Read the definitive tables back from a build directory.

    This is what makes the legacy export a projection: it reads what shipped, never ``chunks/``.
    """
    return {name: read_table(d, name) for name in _TABLE_ORDER}


def write_all(tables: dict[str, pl.DataFrame], out: Path) -> dict[str, Path]:
    """Every new-format member: three tables, the joined view, and the generated schema."""
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    frames = {**{n: tables[n].select(cols) for n, cols in _TABLE_ORDER.items()},
              "vdjdb": joined(tables)}
    for name, frame in frames.items():
        # The one place the tidy names are applied. `strict=False` is wrong here: every column of
        # these tables is declared, so an unmapped one is a registry gap and should raise.
        frame = frame.rename({c: _TIDY[c] for c in frame.columns})
        frame.write_parquet(out / f"{name}.parquet")
        frame.write_csv(out / f"{name}.tsv", separator="\t", line_terminator="\n")
        written[f"{name}.parquet"] = out / f"{name}.parquet"
        written[f"{name}.tsv"] = out / f"{name}.tsv"

    # Dtypes are read off the frames rather than declared, so the schema cannot claim a type the
    # files do not have.
    # Only the three fact tables: the joined view's columns are theirs, and `vdjdb` already names
    # the legacy 22-column table in the registry.
    dtypes = {n: {c: str(t) for c, t in frames[n].schema.items()} for n in _TABLE_ORDER}
    (out / "vdjdb.schema.json").write_text(schema_json(dtypes))
    written["vdjdb.schema.json"] = out / "vdjdb.schema.json"
    return written
