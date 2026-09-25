"""The legacy export: a join and a pivot off the definitive tables.

Nothing here assembles anything. ``records`` and ``chains`` (:mod:`vdjdb.assemble.tables`) are the
database; this module projects them back into the three shapes ``vdjdb-web`` and standalone clients
have always read, and it is the **only** place that knows about them.

That direction matters. The legacy files are untidy in three specific ways -- paired alpha/beta
columns, record fields duplicated per chain, and JSON blobs -- and each is undone here by exactly
one operation:

=========================  ==================================================================
``vdjdb_full.txt``         ``chains`` pivoted **wide** on ``gene``, joined to ``records``
``vdjdb.txt``              ``chains`` joined to ``records``, blobs reassembled
``vdjdb.slim.txt``         ``vdjdb.txt`` grouped, non-key fields collapsed to sorted sets
=========================  ==================================================================

``vdjdb-web`` parses these positionally (CLAUDE.md hard rule 1), so column order is a contract.
Every deliberate quirk reproduced here carries a comment saying why; anything without one is a bug.
"""
from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from ..assemble.master import FIX_FIELDS
from ..schema import FULL_COLUMNS, META_COLUMNS, METHOD_COLUMNS, SLIM_COLUMNS, VDJDB_COLUMNS
from ..schema.fields import render_meta, render_slim_meta

#: Fed to every writer. The release contains no quoted field -- the JSON cells contain ``"``
#: characters, but a default CSV writer would *wrap* them and double the inner quotes, changing
#: every JSON cell in the database.
_CSV = {"separator": "\t", "quote_style": "never", "line_terminator": "\n",
        "include_header": True}

#: Deduplication key for the two ``.found`` counters. Eleven columns: no reference, no donor --
#: the same clonotype seen twice is one sample.
SAMPLE_SIGNATURE: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "d.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope",
)

SLIM_GROUP: tuple[str, ...] = ("gene", "cdr3", "species",
                               "antigen.epitope", "antigen.gene", "antigen.species")

_CHAIN_SUFFIX = {"TRA": "alpha", "TRB": "beta"}


# ---------------------------------------------------------------------------------------------
# chains -> the wide, paired shape
# ---------------------------------------------------------------------------------------------

def widen_chains(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """One row per record, with both chains side by side.

    The inverse of :func:`vdjdb.assemble.tables.build_chains`. Half the cells of the result are
    empty by construction, which is why the stored shape is long.
    """
    out = records
    for tag, suffix in _CHAIN_SUFFIX.items():
        side = chains.filter(pl.col("gene") == tag).select(
            "record_id",
            pl.col("cdr3").alias(f"cdr3.{suffix}"),
            pl.col("v.segm").alias(f"v.{suffix}"),
            pl.col("j.segm").alias(f"j.{suffix}"),
            pl.col("d.segm").alias(f"d.{suffix}"),
            *(pl.col(_CHAIN_COL[key]).alias(f"__{key}.{suffix}") for key, _, _ in FIX_FIELDS),
            pl.col("TCR_hash").alias(f"__hash.{suffix}"),
        )
        out = out.join(side, on="record_id", how="left")
    # A join does not promise an order. `complex.id` is a positional counter, so an unordered
    # result would renumber every clone on every run -- exactly the class of defect CLAUDE.md hard
    # rule 7 exists for. Curation order is the one order that means something here.
    out = out.sort("chunk.file", "chunk.row")
    # An absent chain is an empty string, never a null: one missing marker (CLAUDE.md rule 6).
    fills = [f"cdr3.{s}" for s in _CHAIN_SUFFIX.values()]
    fills += [f"{p}.{s}" for s in _CHAIN_SUFFIX.values() for p in ("v", "j", "d")]
    return out.with_columns(
        *(pl.col(c).fill_null("") for c in fills if c in out.columns),
        pl.coalesce("__hash.alpha", "__hash.beta").fill_null("").alias("TCR_hash"),
    ).drop("__hash.alpha", "__hash.beta")


#: ``cdr3fix`` member -> the ``chains`` column holding it.
_CHAIN_COL = {
    "cdr3": "cdr3", "cdr3_old": "cdr3.original", "fixNeeded": "fix.needed", "good": "fix.good",
    "jCanonical": "j.canonical", "jFixType": "j.fix.type", "jId": "j.segm", "jStart": "j.start",
    "vCanonical": "v.canonical", "vEnd": "v.end", "vFixType": "v.fix.type", "vId": "v.segm",
}


def _cdr3fix_blob(suffix: str, *, as_json: bool) -> pl.Expr:
    """Reassemble the ``cdr3fix`` JSON blob from the flat columns, in the historical key order.

    ``vdjdb.txt`` ships JSON; ``vdjdb_full.txt`` ships a Python ``dict`` **repr** -- single quotes,
    so no standard parser reads it. All 122,930 non-empty cells of the 2026-06-03 release are that
    form. It is a defect, reproduced here and fixed in the definitive tables, where every member is
    already a column.
    """
    render = (lambda d: json.dumps(dict(d))) if as_json else (lambda d: repr(dict(d)))
    fields = {key: pl.col(f"__{key}.{suffix}") for key, _, _ in FIX_FIELDS}
    return (
        pl.when(pl.col(f"cdr3.{suffix}") == "").then(pl.lit(""))
        .otherwise(pl.struct(**fields).map_elements(render, return_dtype=pl.Utf8))
    )


def write_full(records: pl.DataFrame, chains: pl.DataFrame, out: Path) -> Path:
    """``vdjdb_full.txt`` -- the 31 chunk columns plus four derived, one row per record."""
    wide = widen_chains(records, chains).with_columns(
        *(_cdr3fix_blob(s, as_json=False).alias(f"cdr3fix.{s}") for s in _CHAIN_SUFFIX.values()),
    )
    path = out / "vdjdb_full.txt"
    wide.select(FULL_COLUMNS).write_csv(path, **_CSV)
    return path


# ---------------------------------------------------------------------------------------------
# records x chains -> the flat, per-chain shape
# ---------------------------------------------------------------------------------------------

def _web_method() -> pl.Expr:
    """A coarse identification class, for fast filtering in the web front end."""
    m = pl.col("method.identification").str.to_lowercase()
    return (
        pl.when(m.str.contains("sort", literal=True)).then(pl.lit("sort"))
        .when(m.str.contains("culture", literal=True) | m.str.contains("cloning", literal=True)
              | m.str.contains("targets", literal=True)).then(pl.lit("culture"))
        .otherwise(pl.lit("other"))
    )


def _web_method_seq() -> pl.Expr:
    s = pl.col("method.sequencing").str.to_lowercase()
    return (
        pl.when(pl.col("method.singlecell") != "").then(pl.lit("singlecell"))
        .when(s.str.contains("sanger", literal=True)).then(pl.lit("sanger"))
        .when(s.str.contains("-seq", literal=True)).then(pl.lit("amplicon"))
        .otherwise(pl.lit("other"))
    )


def _json_map(prefix: str, columns: tuple[str, ...],
              extra: dict[str, pl.Expr] | None = None) -> pl.Expr:
    """A JSON object column with ``json.dumps`` defaults.

    ``struct.json_encode()`` is **not** byte-compatible: it emits ``{"a":"x"}`` where the release
    has ``{"a": "x"}``, and it does not escape non-ASCII, where ``json.dumps`` writes
    ``M158\\u201366``. Measured at 0.25 s for 200k rows, which is not worth a byte-level fight.
    """
    fields = {c.removeprefix(prefix): pl.col(c) for c in columns}
    if extra:
        fields |= extra
    return pl.struct(**fields).map_elements(lambda d: json.dumps(dict(d)), return_dtype=pl.Utf8)


def build_default(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """``vdjdb.txt`` -- one row per chain, paired chains sharing a ``complex.id``.

    ``complex.id`` is assigned **here** and nowhere else: it is a positional counter with no
    meaning outside this file, and the definitive tables express pairing by ``record_id``.
    """
    wide = widen_chains(records, chains)

    # A chain with a sequence but no V or J cannot be placed on a germline, so the record is
    # dropped whole -- both chains -- exactly as the legacy did.
    ok = pl.all_horizontal(
        *[((pl.col(f"cdr3.{g}") == "") | ((pl.col(f"v.{g}") != "") & (pl.col(f"j.{g}") != "")))
          for g in _CHAIN_SUFFIX.values()]
    )
    wide = wide.filter(ok)

    paired = (pl.col("cdr3.alpha") != "") & (pl.col("cdr3.beta") != "")
    wide = wide.with_columns(
        pl.when(paired).then(paired.cum_sum()).otherwise(0).cast(pl.Int64).alias("complex.id"),
        pl.len().over(SAMPLE_SIGNATURE).alias("__samples"),
        pl.col("reference.id").n_unique().over(SAMPLE_SIGNATURE).alias("__studies"),
        _web_method().alias("web.method"),
        _web_method_seq().alias("web.method.seq"),
        *(_cdr3fix_blob(s, as_json=True).alias(f"__fixjson.{s}") for s in _CHAIN_SUFFIX.values()),
    )

    parts = []
    for tag, suffix in _CHAIN_SUFFIX.items():
        present = pl.col(f"cdr3.{suffix}") != ""
        parts.append(
            wide.filter(present).with_columns(
                pl.lit(tag).alias("gene"),
                pl.col(f"cdr3.{suffix}").alias("cdr3"),
                pl.col(f"v.{suffix}").alias("v.segm"),
                pl.col(f"j.{suffix}").alias("j.segm"),
                pl.col(f"__fixjson.{suffix}").alias("cdr3fix"),
                # Canonical anchors: is this CDR3 bounded by the expected Cys and Phe/Trp?
                pl.when(pl.col(f"__vCanonical.{suffix}") & pl.col(f"__jCanonical.{suffix}"))
                .then(pl.lit("no")).otherwise(pl.lit("yes")).alias("web.cdr3fix.nc"),
                # Unmapped: could the CDR3 be placed on a V *and* a J germline? The shipped build
                # tests `cdr3fix["jStart"]` for truthiness where the Groovy tested `jStart > -1`,
                # and -1 is truthy in Python -- so 7,973 of 284,546 unmapped records are labelled
                # mapped. Fixed here; declared as `web-unmp-jstart-minus1` in the ledger.
                pl.when((pl.col(f"__vEnd.{suffix}") > -1) & (pl.col(f"__jStart.{suffix}") > -1))
                .then(pl.lit("no")).otherwise(pl.lit("yes")).alias("web.cdr3fix.unmp"),
            )
        )
    melted = pl.concat(parts, how="vertical")

    return melted.with_columns(
        _json_map("method.", METHOD_COLUMNS).alias("method"),
        _json_map("meta.", META_COLUMNS, {
            "samples.found": pl.col("__samples").cast(pl.Int64),
            "studies.found": pl.col("__studies").cast(pl.Int64),
        }).alias("meta"),
    ).select(VDJDB_COLUMNS)


def build_slim(default_db: pl.DataFrame) -> pl.DataFrame:
    """``vdjdb.slim.txt`` -- one row per (chain, epitope), every other field collapsed.

    Non-key fields become a sorted comma-separated set, which is how a consumer sees "this CDR3 was
    reported against this epitope in these studies, with these V genes". ``vdjdb.score`` is the
    **maximum** rather than a set: a score is a confidence, and the best evidence wins.
    """
    fix = pl.col("cdr3fix").str.json_decode(
        dtype=pl.Struct([pl.Field("vEnd", pl.Int64), pl.Field("jStart", pl.Int64)]))
    df = default_db.with_columns(
        pl.when(pl.col("cdr3fix") == "").then(pl.lit(""))
        .otherwise(fix.struct.field("jStart").cast(pl.Utf8)).alias("j.start"),
        pl.when(pl.col("cdr3fix") == "").then(pl.lit(""))
        .otherwise(fix.struct.field("vEnd").cast(pl.Utf8)).alias("v.end"),
    )
    collapse = [c for c in SLIM_COLUMNS if c not in SLIM_GROUP and c != "vdjdb.score"]
    return (
        df.group_by(SLIM_GROUP)
        .agg(pl.col("vdjdb.score").max(),
             *[pl.col(c).cast(pl.Utf8).unique().sort().str.join(",").alias(c) for c in collapse])
        .sort(SLIM_GROUP)          # pandas `groupby` sorts by key; reproduce it explicitly
        .select(SLIM_COLUMNS)
    )


def write_meta(out: Path) -> list[Path]:
    """The two metadata files, **generated** rather than tracked, so they cannot drift again."""
    a, b = out / "vdjdb.meta.txt", out / "vdjdb.slim.meta.txt"
    a.write_text(render_meta("vdjdb"))
    b.write_text(render_slim_meta("slim"))
    return [a, b]


def write_all(tables: dict[str, pl.DataFrame], out: Path) -> dict[str, Path]:
    """Every legacy member, from the definitive tables."""
    out.mkdir(parents=True, exist_ok=True)
    records, chains = tables["records"], tables["chains"]
    default_db = build_default(records, chains)
    slim = build_slim(default_db)
    paths = {"vdjdb_full.txt": write_full(records, chains, out)}
    default_db.write_csv(out / "vdjdb.txt", **_CSV)
    slim.write_csv(out / "vdjdb.slim.txt", **_CSV)
    paths["vdjdb.txt"] = out / "vdjdb.txt"
    paths["vdjdb.slim.txt"] = out / "vdjdb.slim.txt"
    for p in write_meta(out):
        paths[p.name] = p
    return paths
