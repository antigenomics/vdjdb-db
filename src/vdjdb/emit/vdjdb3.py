"""The new VDJdb format: the definitive tables as they ship, plus the joined convenience view.

Three fact tables and one view, parquet with a TSV projection of each (`docs/outputs.md` §3):

=====================  =========================================================================
``records``            one row per curated record, PK ``record_id``
``chains``             one row per TCR chain, PK ``(record_id, gene)``
``evidence``           one row per piece of support, PK ``(record_id, evidence_id)``
``vdjdb``              ``records`` |><| ``chains`` |><| ``evidence``, the denormalised view
=====================  =========================================================================

The first three are written, not built: :mod:`vdjdb.assemble.tables` produced them and this module
only chooses a file format. ``vdjdb`` is the one thing derived here, and it is derived every time --
a consumer who edits it is editing a cache, which is exactly the property the legacy format lacked.

Unlike :mod:`vdjdb.emit.legacy`, nothing here reproduces a byte layout, so the writers are plain:
default quoting rather than ``quote_style="never"`` (a field containing a tab must survive, not
corrupt the row), and no JSON blobs at all.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..schema import CHAIN_COLUMNS, EVIDENCE_TABLE_COLUMNS, RECORD_COLUMNS, schema_json

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
#: without a producer yet is ``false`` -- an honest "no evidence of this kind", not a missing column.
#: ``same.study`` has no producer at all: method-level self-validation is phase 9's.
VIEW_EVIDENCE_COLUMNS: tuple[str, ...] = (*EVIDENCE_VIEW.values(), "evidence.validation.same.study")

_TABLE_ORDER: dict[str, tuple[str, ...]] = {
    "records": RECORD_COLUMNS, "chains": CHAIN_COLUMNS, "evidence": EVIDENCE_TABLE_COLUMNS,
}


def joined(tables: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """``records`` joined to ``chains``, with each evidence type pivoted to a boolean column.

    One row per chain -- the level ``vdjdb.txt`` is written at -- so this is the table a user who
    just wants "one flat thing" should read, and it is the one table here that is purely derived.
    """
    records, chains, evidence = tables["records"], tables["chains"], tables["evidence"]
    out = chains.join(records, on="record_id", how="left")

    ev = evidence.with_columns(
        pl.col("evidence_type").replace_strict(EVIDENCE_VIEW, default=None).alias("__col")
    ).drop_nulls("__col")
    if ev.select(pl.col("gene").fill_null("").eq("").any()).item():
        # Record-level evidence (a structure covers the receptor, not one chain) needs a second
        # join on `record_id` alone. Nothing produces it yet; refuse rather than drop it silently.
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


def read_tables(d: Path) -> dict[str, pl.DataFrame]:
    """Read the definitive tables back from a build directory.

    This is what makes the legacy export a *projection*: it reads what shipped, never ``chunks/``.
    """
    return {name: pl.read_parquet(d / f"{name}.parquet") for name in _TABLE_ORDER}


def write_all(tables: dict[str, pl.DataFrame], out: Path) -> dict[str, Path]:
    """Every new-format member: three tables, the joined view, and the generated schema."""
    out.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    frames = {**{n: tables[n].select(cols) for n, cols in _TABLE_ORDER.items()},
              "vdjdb": joined(tables)}
    for name, frame in frames.items():
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
