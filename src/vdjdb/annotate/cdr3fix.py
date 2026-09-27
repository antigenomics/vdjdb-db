"""CDR3 markup and repair, via ``arda.cdr3fix``.

Replaces the k-mer scanner vendored from the Groovy-era build
(``annotate/_legacy_fixer/``). ``arda`` aligns the germline-templated run with a semi-global
Needleman-Wunsch anchored at [FW]118 with free end gaps, where the legacy used a longest-hit k-mer
scan, so it extends the match further into the junction and reports it starting earlier.

Measured on 20,000 rows of the 2026-06-03 release (``random.seed(42)``), agreement with the shipped
``cdr3fix``:

=================  ===========  ==============================================================
``cdr3``           97.81 %
``vFixType``       97.95 %
``vEnd``           96.53 %
``good``           95.68 %
``jFixType``       95.44 %
``jStart``         91.02 %      of the 1,797 disagreements: 556 VDJdb-unmapped that arda maps,
                                1,241 where both map and arda is smaller in every case
                                (mode -2, range -1..-8), and 0 coverage regressions
=================  ===========  ==============================================================

The deviation is one-directional: arda never loses a mapping the legacy had.

``Cdr3Markup.to_cdr3fix()`` emits VDJdb's JSON keys and fix-type names verbatim, so the flat
columns below hold the same values the legacy ones did.

Species. VDJdb spells species ``HomoSapiens``; arda wants ``human``. ``arda.cdr3fix``'s own
``VDJDB_SPECIES`` map is the translation and the authority -- do not hand-roll one.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from ..assemble.master import _FIX_DTYPES, FIX_FIELDS


def ensure_reference() -> Path:
    """Make sure arda resolves a germline reference, and raise if it cannot.

    arda decides it is running from a source checkout by walking up from its own ``__file__``
    looking for a directory with ``database/`` and a project marker. Installed into this project's
    ``.venv``, that walk reaches this repository, which has both -- so arda points at
    ``vdjdb-db/database/vdj``, which does not exist. Before arda 2.28 ``load_anchors`` then returned
    an empty dict rather than raising, and all 191,447 CDR3s came back ``FailedBadSegment`` with
    ``vEnd = -1``: a misconfiguration becomes a database of wrong annotations, with no error.

    Fixed upstream (arda ``_source_root`` now requires ``database/vdj``). Until that release is
    pinned, point ``$ARDA_HOME`` at the per-user cache, fetching the reference into it first if it
    is not there, because arda's own auto-fetch fires only when ``_source_root()`` is ``None`` and
    here it is not. Without that, a clean machine has no reference at all: the first CI run of
    ``build.yml`` failed this way.
    """
    from arda.cdr3fix import load_anchors
    from arda.paths import cache_root, database_dir

    if (database_dir() / "vdj").is_dir() and load_anchors("human"):
        return database_dir()

    cache = cache_root()
    if not (cache / "database" / "vdj").is_dir() and "ARDA_NO_AUTO_FETCH" not in os.environ:
        # arda auto-fetches the reference, but ONLY when `_source_root()` is None -- and here it
        # is not, because the walk finds this repository. So the fetch that would have happened on
        # a plain install never fires, and a clean machine (a CI runner, a new checkout) has no
        # reference at all. Fetching it here lets the build run from nothing.
        #
        # This is a download of an input, not a cache of a result: hard rule 9's explicit
        # carve-out. The germline reference is data arriving, not a result being remembered.
        from arda._database_fetch import fetch_database

        fetch_database(cache / "database")

    if (cache / "database" / "vdj").is_dir():
        os.environ["ARDA_HOME"] = str(cache)
        for fn in (load_anchors, __import__("arda.paths", fromlist=["x"])._source_root,
                   __import__("arda.paths", fromlist=["x"]).database_dir):
            fn.cache_clear()
        if load_anchors("human"):
            return database_dir()

    raise RuntimeError(
        f"arda has no germline reference: database_dir() is {database_dir()} and it has no "
        f"usable vdj/. Set $ARDA_HOME to an arda checkout or to a cache root containing "
        f"database/vdj, or unset $ARDA_NO_AUTO_FETCH so arda can fetch it."
    )

#: ``(species, cdr3, v, j)`` -- one markup per distinct key, never per row. Measured at 27 us per
#: record, so the 191,447 distinct keys cost ~6 s; the join back is free.
KEY: tuple[str, ...] = ("species", "cdr3", "v", "j")

#: Residues arda may substitute to make a CDR3 conform to the germline it was handed. Zero.
#:
#: The legacy default was 1, and arda aligns more sensitively, so keeping it would have rewritten
#: 4,486 curated human CDR3s rather than 649 -- and the rewrites are wrong for this database.
#: `CAAADSWGKLQF` with `TRAJ24*01` becomes `CAAADSWGKLEF`, because `WGKLEF` is what `*01` encodes.
#: But `WGKLQF` is the `*02` signature: the sequence is right and the allele call is wrong
#: (issue #327, where 66 % of explicit `*01` calls show the `*02` motif). Substituting the residue
#: destroys the evidence that would fix the call.
#:
#: Measured, it costs nothing: both ends map on 165,223 of 174,630 distinct human keys at
#: either setting, and `good` is marginally higher at 0 (164,925 against 164,899). Trimming and
#: extending are unaffected -- 2,588 against 2,620 -- because those repair a truncated sequence
#: rather than contradicting a reported one.
MAX_REPLACE = 0


def markup(keys: pl.DataFrame, gene: str | None = None) -> pl.DataFrame:
    """Mark up a frame of distinct ``(species, cdr3, v, j)`` keys.

    One :func:`arda.cdr3fix.markup_records` call per organism -- anchors are loaded and cached once
    per organism, so grouping the call that way is the difference between one load and 100,000.
    """
    from arda.cdr3fix import VDJDB_SPECIES, markup_records

    ensure_reference()
    # The guessed segments go in their own columns: `v` and `j` are the join key back to the
    # table, so overwriting them would make the lookup miss -- measured, it fanned vdjdb_full.txt
    # out by 3,266 rows.
    filled = guess_missing_segments(keys, gene)
    out: list[pl.DataFrame] = []
    for (species,) in filled.select("species").unique().sort("species").iter_rows():
        organism = VDJDB_SPECIES.get(species.lower())
        part = filled.filter(pl.col("species") == species)
        if organism is None:
            # An unknown species is not a reason to drop records: mark them unmapped and let the
            # QC rules complain about the species, which is the actual defect.
            out.append(part.with_columns(_unmapped(part)))
            continue
        records = markup_records(part, v="__gv", j="__gj", organism=organism,
                                 max_replace=MAX_REPLACE)
        fixes = [r.to_cdr3fix() for r in records]
        out.append(part.with_columns(
            *(pl.Series(tmp, [f[key] for f in fixes], dtype=_FIX_DTYPES[ty])
              for key, tmp, ty in FIX_FIELDS)
        ).with_columns(
            # arda returns an empty segment id when it cannot resolve the call; the legacy kept
            # the closest match it had. Dropping the call as well as the coordinates would fail
            # the legacy "a CDR3 needs a V and a J" filter and cost 11,619 rows of vdjdb.txt --
            # a coverage regression, which phase 5 must not produce. The coordinates stay -1,
            # recording that nothing was located.
            pl.when(pl.col("__v") == "").then(pl.col("__gv")).otherwise(pl.col("__v")).alias("__v"),
            pl.when(pl.col("__j") == "").then(pl.col("__gj")).otherwise(pl.col("__j")).alias("__j"),
        ))
    return pl.concat(out, how="vertical").drop("__gv", "__gj").sort(KEY)


def guess_missing_segments(keys: pl.DataFrame, gene: str | None = None) -> pl.DataFrame:
    """Name a V or J for rows that have none, using the legacy k-mer guesser.

    arda repairs a CDR3; it does not guess a segment. Given a blank ``v``, ``markup_records``
    reports ``FailedBadSegment`` rather than proposing one, and the record then fails the legacy
    build's "a CDR3 needs a V and a J" filter. Measured: swapping both halves at once dropped
    13,844 of 284,546 rows from ``vdjdb.txt`` -- a coverage regression, which the phase 5
    acceptance criterion forbids.

    So the guesser stays on the legacy k-mer scan for now. Replacing it is its own deviation, with
    its own measurement: issue #462 and ROADMAP phase 8, where candidate segments are scored by
    OLGA Pgen instead of longest-hit.
    """
    blank = (pl.col("v") == "") | (pl.col("j") == "")
    if not keys.filter(blank).height:
        return keys.with_columns(pl.col("v").alias("__gv"), pl.col("j").alias("__gj"))

    from .._legacy_guess import guess_segments

    return guess_segments(keys, gene)


def _unmapped(part: pl.DataFrame) -> list[pl.Expr]:
    """What an unmappable record looks like: the sequence unchanged, nothing located."""
    blank = {"cdr3": pl.col("cdr3"), "cdr3_old": pl.col("cdr3"), "fixNeeded": pl.lit(False),
             "good": pl.lit(False), "jCanonical": pl.lit(False),
             "jFixType": pl.lit("FailedBadSegment"), "jId": pl.lit(""), "jStart": pl.lit(-1),
             "vCanonical": pl.lit(False), "vEnd": pl.lit(-1),
             "vFixType": pl.lit("FailedBadSegment"), "vId": pl.lit("")}
    return [blank[key].cast(_FIX_DTYPES[ty]).alias(tmp) for key, tmp, ty in FIX_FIELDS]
