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

    from .segments import propose

    ensure_reference()
    # The proposal goes in its own columns: `v` and `j` are the join key back to the table, so
    # overwriting them would make the lookup miss -- measured, it fanned vdjdb_full.txt out by 3,266
    # rows. `__mj` is what the repair runs against and `__gj` is what the table reports; they differ
    # only in being blank where the record named its own segment.
    filled = propose(keys, gene).with_columns(
        # **Only the J proposal reaches the repair, and through it the shipped call.** The two sides
        # are not equally knowable from a junction and the measurement says so: a J is recovered at
        # 93.6-97.5 % gene-level accuracy, because its germline templates a distinctive 3' motif,
        # while a V is recovered at 23.8-50.1 %, because TRBV contributes only a few residues and
        # most of them template the same `CAS`. Feeding the V proposal in was tried: all 706 chains
        # it reached came back `NoFixNeeded` with no repair proposed, so arda confirmed that *a*
        # germline fits without discriminating between the many that fit equally, and `v.segm` would
        # then carry one of them as if the publication had reported it.
        #
        # So the V proposal is reported as `v.inferred` and nothing else, `v.segm` stays blank where
        # the curator left it blank and arda resolved nothing - which is what every release has
        # shipped - and the V boundary question is answered by `v.end.inferred`.
        pl.col("v").alias("__mv"), pl.col("__gj").alias("__mj"))
    out: list[pl.DataFrame] = []
    for (species,) in filled.select("species").unique().sort("species").iter_rows():
        organism = VDJDB_SPECIES.get(species.lower())
        part = filled.filter(pl.col("species") == species)
        if organism is None:
            # An unknown species is not a reason to drop records: mark them unmapped and let the
            # QC rules complain about the species, which is the actual defect.
            out.append(part.with_columns(_unmapped(part)).with_columns(
                # Same schema as the mapped branch: arda named nothing, so it proposes nothing.
                pl.lit("").alias("__varda"), pl.lit("").alias("__jarda"),
                pl.col("v").alias("__v"), pl.col("j").alias("__j"),
            ))
            continue
        records = markup_records(part, v="__mv", j="__mj", organism=organism,
                                 max_replace=MAX_REPLACE)
        fixes = [r.to_cdr3fix() for r in records]
        out.append(part.with_columns(
            *(pl.Series(tmp, [f[key] for f in fixes], dtype=_FIX_DTYPES[ty])
              for key, tmp, ty in FIX_FIELDS)
        ).with_columns(
            # arda's own call is kept as evidence, under its own name, and is what `chains` reports
            # as `v.segm.arda` / `j.segm.arda`. It is not what ships: see the two rules below.
            pl.col("__v").str.replace_all(";", ",").alias("__varda"),
            pl.col("__j").str.replace_all(";", ",").alias("__jarda"),
        ).with_columns(
            # **The engine's allele-resolved call ships, and that is not new.** The markup engine
            # names the allele it aligned against, so a bare `TRBV12-3` comes back `TRBV12-3*01`.
            # Measured before assuming it: the 2026-06-03 release carries an allele suffix on
            # **282,277 of 284,546** `v.segm` cells (99.2 %) and 282,763 `j.segm` (99.4 %), while
            # only **37.0 %** of submitted `v.beta` and 40.2 % of `j.beta` carry one. So every
            # release VDJdb has shipped already published the fixer's resolved allele rather than
            # the curator's bare call, and keeping the bare call instead was tried and moved
            # 183,345 of 284,546 rows -- a far bigger change than the swap it was meant to avoid.
            # `chains` records what was submitted next to it, so neither is lost.
            #
            # Where the engine cannot resolve the call the guesser's stands: dropping it would fail
            # the legacy "a CDR3 needs a V and a J" filter and cost 11,619 rows of `vdjdb.txt`. The
            # coordinates then stay -1, recording that nothing was located.
            # `__mv` / `__mj`, not `__gv` / `__gj`: the fallback is the call the repair was given,
            # which on the V side is the record's own and on the J side includes the proposal.
            pl.when(pl.col("__varda") != "").then(pl.col("__varda"))
              .otherwise(pl.col("__mv")).alias("__v"),
            pl.when(pl.col("__jarda") != "").then(pl.col("__jarda"))
              .otherwise(pl.col("__mj")).alias("__j"),
        ))
    return (pl.concat(out, how="vertical")
            .with_columns(
                # Reported as `v.inferred` / `j.inferred`: the proposal alone, blank wherever the
                # record named the segment itself, so a proposal never sits beside a curated call.
                *(pl.when(pl.col(call) == "").then(pl.col(tmp)).otherwise(pl.lit("")).alias(tmp)
                  for call, tmp in (("v", "__gv"), ("j", "__gj"))))
            .drop("__mv", "__mj").sort(KEY))


def _unmapped(part: pl.DataFrame) -> list[pl.Expr]:
    """What an unmappable record looks like: the sequence unchanged, nothing located."""
    blank = {"cdr3": pl.col("cdr3"), "cdr3_old": pl.col("cdr3"), "fixNeeded": pl.lit(False),
             "good": pl.lit(False), "jCanonical": pl.lit(False),
             "jFixType": pl.lit("FailedBadSegment"), "jId": pl.lit(""), "jStart": pl.lit(-1),
             "vCanonical": pl.lit(False), "vEnd": pl.lit(-1),
             "vFixType": pl.lit("FailedBadSegment"), "vId": pl.lit("")}
    return [blank[key].cast(_FIX_DTYPES[ty]).alias(tmp) for key, tmp, ty in FIX_FIELDS]
