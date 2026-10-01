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


#: VDJdb's chain word -> the locus a proposal must resolve to before it may ship. arda 2.36 resolves
#: the locus itself from the junction, which is what lets a record naming neither V nor J get a call
#: at all - but this schema has an alpha column and a beta column and nothing else, so a proposal
#: that lands on a third locus is evidence about the record rather than a call for it. Measured: of
#: 461 keys naming neither side, arda's locus agrees with the column the record was filed under on
#: **457 (99.13 %)**, and all four disagreements are `CACD...DKLIF` - TRDV2's own anchor and TRDJ1's
#: own ending, filed as alpha because there is nowhere else to put it.
LOCI: dict[str, str] = {"alpha": "TRA", "beta": "TRB"}


def markup(keys: pl.DataFrame, gene: str | None = None, *, recall_j: bool = True) -> pl.DataFrame:
    """Mark up a frame of distinct ``(species, cdr3, v, j)`` keys.

    One :func:`arda.cdr3fix.markup_records` call for the whole frame: it reads ``species`` per row
    and loads each organism's anchors once, so there is nothing to group by here.

    **arda proposes a call for a side the submission left blank, and this used to be our job.**
    ``annotate/segments.py`` held a k-mer candidate lookup with a Pgen tiebreak for exactly that,
    and arda 2.34 took over the one-sided cases while 2.36 added the locus, so all 3,130 blank-call
    keys are now answered by the engine that repairs the junction rather than by a second one
    beside it. That module is deleted: the germline a repair ran against and the call reported next
    to it could disagree while they were two computations, and now they cannot.
    """
    from arda.cdr3fix import markup_records

    ensure_reference()
    out = _frame(keys, markup_records(keys, max_replace=MAX_REPLACE), gene)
    if recall_j:
        out = _recall_j(out, gene)
    return out.sort(KEY)


def _frame(keys: pl.DataFrame, records: list, gene: str | None) -> pl.DataFrame:
    """The shipped columns for ``keys``, from one ``Cdr3Markup`` per key and in the same order."""
    fixes = [r.to_cdr3fix() for r in records]
    locus = LOCI.get(gene or "")
    # A proposal is accepted only for the locus this column can hold; see LOCI.
    usable = [locus is None or r.locus == locus for r in records]
    out = keys.with_columns(
        *(pl.Series(tmp, [f[key] for f in fixes], dtype=_FIX_DTYPES[ty])
          for key, tmp, ty in FIX_FIELDS),
        pl.Series("__pv", [("V" in r.proposed) and u for r, u in zip(records, usable, strict=True)],
                  dtype=pl.Boolean),
        pl.Series("__pj", [("J" in r.proposed) and u for r, u in zip(records, usable, strict=True)],
                  dtype=pl.Boolean),
        pl.Series("__locus", [r.locus or "" for r in records], dtype=pl.Utf8),
    ).with_columns(
        # arda's own call is kept as evidence, under its own name, and is what `chains` reports as
        # `v.segm.arda` / `j.segm.arda`. It is not what ships: see the two rules below.
        pl.col("__v").str.replace_all(";", ",").alias("__varda"),
        pl.col("__j").str.replace_all(";", ",").alias("__jarda"),
    ).with_columns(
        # **Only the J proposal reaches the shipped call.** The two sides are not equally knowable
        # from a junction and the measurement says so: a J is recovered at 93.6-97.5 % gene-level
        # accuracy, because its germline templates a distinctive 3' motif, while a V is recovered at
        # 23.8-50.1 %, because TRBV contributes only a few residues and most of them template the
        # same `CAS`. So a proposed V is reported as `v.inferred` and nowhere else, `v.segm` stays
        # blank where the curator left it blank - which is what every release has shipped - and the
        # V boundary question is answered by `v.end.inferred`.
        #
        # **The engine's allele-resolved call ships where the record named the segment, and that is
        # not new.** The markup engine names the allele it aligned against, so a bare `TRBV12-3`
        # comes back `TRBV12-3*01`. Measured before assuming it: the 2026-06-03 release carries an
        # allele suffix on **282,277 of 284,546** `v.segm` cells (99.2 %) and 282,763 `j.segm`
        # (99.4 %), while only **37.0 %** of submitted `v.beta` and 40.2 % of `j.beta` carry one. So
        # every release VDJdb has shipped already published the fixer's resolved allele rather than
        # the curator's bare call. `chains` records what was submitted next to it, so neither is
        # lost. Where the engine resolves nothing the submitted call stands: dropping it would fail
        # the legacy "a CDR3 needs a V and a J" filter and cost 11,619 rows of `vdjdb.txt`, and the
        # coordinates then stay -1, recording that nothing was located.
        _shipped_call("v", "__varda", proposed="__pv").alias("__v"),
        _shipped_call("j", "__jarda", proposed="__pj").alias("__j"),
        # Reported as `v.inferred` / `j.inferred`: the proposal alone, blank wherever the record
        # named the segment itself, so a proposal never sits beside a curated call.
        pl.when(pl.col("__pv")).then(pl.col("__varda")).otherwise(pl.lit("")).alias("__gv"),
        pl.when(pl.col("__pj")).then(pl.col("__jarda")).otherwise(pl.lit("")).alias("__gj"),
    )
    return out.drop("__pv", "__pj", "__locus")


def _recall_j(out: pl.DataFrame, gene: str | None) -> pl.DataFrame:
    """Replace a J call the junction does not support by the one gene that does (#681).

    **The submitted call is untouched and the shipped one is fixed.** `j.segm.submitted` carries
    what the paper reported and `j.segm.arda` what the engine called; this changes `j.segm` and the
    metadata about it - `jId`, `jStart`, `jCanonical` - and nothing about the junction, which the rule
    never alters. The rule is :func:`vdjdb.curate.jcalls.recall` and it runs on the call that would
    ship, so the engine's own re-call stays wherever it already satisfies the rule.

    `jFixType` stays what the engine reported and `good` stays as it was: the rule re-calls a gene,
    it does not repair a residue.

    The engine cannot be asked to mark the junction up again under the rule's gene. `markup_cdr3`
    answers with the J its own gapped alignment prefers, which differs from the rule's gene on 443
    shipped chains, and an anchor table restricted to one gene raises inside arda. So the placement
    is :func:`vdjdb.curate.jcalls.place`, for the ~1,000 keys the rule changes.
    """
    from ..curate.jcalls import place, recall

    locus = LOCI.get(gene or "")
    if locus is None or out.is_empty():
        return out
    pairs = out.filter(pl.col("__j") != "").select("species", "__cdr3", "__j").unique()
    rows = []
    for sp, cdr3, j in pairs.iter_rows():
        new = recall(sp, locus, cdr3, j)
        if new is not None:
            allele, start, canonical = place(sp, cdr3, new)
            rows.append((sp, cdr3, j, allele, start, canonical))
    if not rows:
        return out
    table = pl.DataFrame(rows, schema={"species": pl.Utf8, "__cdr3": pl.Utf8, "__j": pl.Utf8,
                                       "__jn": pl.Utf8, "__jsn": pl.Int64, "__jcn": pl.Boolean},
                         orient="row")
    hit = pl.col("__jn").is_not_null()
    return (out.join(table, on=["species", "__cdr3", "__j"], how="left")
               .with_columns(pl.when(hit).then(pl.col("__jn")).otherwise(pl.col("__j")).alias("__j"),
                             pl.when(hit).then(pl.col("__jsn")).otherwise(pl.col("__jstart"))
                               .alias("__jstart"),
                             pl.when(hit).then(pl.col("__jcn")).otherwise(pl.col("__jcanon"))
                               .alias("__jcanon"))
               .drop("__jn", "__jsn", "__jcn"))


def _shipped_call(submitted: str, engine: str, *, proposed: str) -> pl.Expr:
    """Which of the submitted and the engine's call ships, per row.

    **The engine's call ships wherever it has one, and that is not new.** It names the allele it
    aligned against, so a bare `TRBV12-3` comes back `TRBV12-3*01` and a family name `TRBV20` comes
    back `TRBV20-1*01`; every release VDJdb has shipped has published those. The 2026-06-03 release
    carries an allele suffix on **282,277 of 284,546** `v.segm` cells (99.2 %) and 282,763 `j.segm`
    (99.4 %), against **37.0 %** of submitted `v.beta` and 40.2 % of `j.beta`. Where the engine
    resolves nothing the submitted call stands: dropping it would fail the legacy "a CDR3 needs a V
    and a J" filter and cost 11,619 rows of `vdjdb.txt`, and the coordinates then stay -1.

    arda 2.36 also **re-calls the gene** where the junction contradicts the submission - `TRBV10-3`
    on a `CASS...` junction comes back `TRBV19*01`, because `TRBV10-3` templates `CAIS` and the
    sequence is the evidence. Restricting the engine to allele-and-family *resolution* was built and
    measured, and it is worse on every axis, so it is not what ships:

    =======================================  ==============  ============  =============  ==========
    shipped call                             J gene correct  no allele     `j-calls.tsv`  row delta
    =======================================  ==============  ============  =============  ==========
    the engine's, wherever it has one        **99.27 %**     **67**        **273**        symmetric
    only where it resolves the submission    98.67 %         129           688            2,601 lost
    =======================================  ==============  ============  =============  ==========

    J gene against `isalgo/airr_control`'s nucleotide-established calls over 4,372 human TRB keys;
    "no allele" is shipped cells with no `*` suffix; `j-calls.tsv` is the count of J calls their own
    junction contradicts. Gating it left *more* under-specified calls shipping and a bigger
    curation queue, and made the row buckets asymmetric - rows lost rather than moved. What the
    curator reported is kept beside it as `v.segm.submitted` / `j.segm.submitted`, and arda's own as
    `v.segm.arda` / `j.segm.arda`, so a reader can tell a markup decision from a curation one.

    **A proposed V never ships**, and that is the one gate. The two sides are not equally knowable
    from a junction: a J is recovered at 93.6-97.5 % gene-level accuracy, because its germline
    templates a distinctive 3' motif, while a V is recovered at 23.8-50.1 %, because TRBV
    contributes only a few residues and most of them template the same `CAS`. So a proposed V is
    reported as `v.inferred` and nowhere else, and `v.segm` stays blank where the curator left it
    blank - which is what every release has shipped.
    """
    sub, eng = pl.col(submitted), pl.col(engine)
    return (pl.when(pl.col(proposed) & (submitted == "v")).then(pl.lit(""))
              .when(eng != "").then(eng).otherwise(sub))


def _unmapped(part: pl.DataFrame) -> list[pl.Expr]:
    """What an unmappable record looks like: the sequence unchanged, nothing located."""
    blank = {"cdr3": pl.col("cdr3"), "cdr3_old": pl.col("cdr3"), "fixNeeded": pl.lit(False),
             "good": pl.lit(False), "jCanonical": pl.lit(False),
             "jFixType": pl.lit("FailedBadSegment"), "jId": pl.lit(""), "jStart": pl.lit(-1),
             "vCanonical": pl.lit(False), "vEnd": pl.lit(-1),
             "vFixType": pl.lit("FailedBadSegment"), "vId": pl.lit("")}
    return [blank[key].cast(_FIX_DTYPES[ty]).alias(tmp) for key, tmp, ty in FIX_FIELDS]
