"""The VDJdb confidence score, as polars expressions.

Ported from ``py_src/ScoreFactory.py``, which built its score map with ``iterrows()`` over ~192k
rows and then evaluated it again per row with ``master_table.T.apply``. Every rule below is a
column-level expression instead, and the per-signature maximum is one window.

The score answers "how much should a reader trust this specificity annotation", 0-3:

* **3** if a solved structure is attached -- nothing else can outrank direct evidence;
* otherwise ``min(sequencing confidence, specificity confidence)``, so a beautifully sequenced
  clonotype with a weak specificity assay scores as low as the assay, and vice versa.

Three behaviours are preserved exactly because the released scores depend on them:

* ``method.frequency`` is parsed as ``n/m``, ``x%`` or a bare float, and ``n`` alone is the cell
  count -- ``2/47`` means two cells of forty-seven, and a ``%`` form therefore has no cell count;
* the score is a maximum over the 11-column sample signature, not a per-row value. The same
  clonotype assayed twice takes the better of the two, which is why the score cannot be computed
  before the full table is assembled;
* every method term matches a **substring** of the method field, never a token -- see
  :func:`_contains_any`, which records what that costs and why it stays.

The structure term at the top of this docstring is where the score's one live defect sits.
``meta.structure.id`` is read as "a solved structure is attached" by testing that the field is
non-empty, and 2,765 chunk rows hold a figure or table reference there rather than a PDB id.
Measured 2026-09-28 (vdjdb-db#402): **4,817 of the 6,007 records that score 3 owe it to one of
those values** - 3,752 would fall to 0 and 1,065 to 1 - so four of every five top-confidence
records rest on a reference a reader cannot check. The guard is the ``structure id is not a PDB
id`` QC rule; the repair is a curation decision, so nothing here changes.
"""
from __future__ import annotations

import polars as pl

#: The score signature. Not ``CHUNK_DEDUP_KEY``: it has no reference or donor fields, so the same
#: clonotype seen in two studies shares one score. Both keys are called ``SIGNATURE_COLS`` in the
#: legacy code, which is a source of confusion.
SCORE_SIGNATURE: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope",
)

_SORT_BASED = ("sort", "beads", "separation", "stain")
_STIMULATION_BASED = ("targets",)
_CULTURE_BASED = ("culture", "cloning")


def _lower(col: str) -> pl.Expr:
    return pl.col(col).str.strip_chars().str.to_lowercase()


def _contains_any(col: str, needles: tuple[str, ...]) -> pl.Expr:
    """Substring, not token. That is the released semantics and it must not be tightened here.

    Every method term in this module matches a substring of a comma-joined method field, because
    the legacy scorer did and the released scores encode it. It is the right shape for
    `antigen-loaded-targets` matching `targets`, and the wrong shape for a word that contains
    another word: `indirect` contains `direct`, so a value spelled that way would take the
    `_high_specificity` ceiling of 3 and force the sequencing term to 3 with it, awarding the
    maximum for the opposite of what it says. That `direct` test is written inline in
    `_high_specificity` rather than routed through here, but it is the same matching rule.

    Measured 2026-09-28 across all 203,348 chunk rows: `method.verification` holds `direct` (338),
    `direct,antigen-loaded-targets` (82) and `antigen-loaded-targets,direct` (2) and nothing else,
    so no value is currently matched loosely and no score is affected. Recorded rather than fixed
    because tightening one needle and not the other six would be inconsistent, and tightening all
    seven diverges from the release the comparison is measured against. If it ever has to change,
    it changes for all of them at once, with a declared rule and a re-measured count.
    """
    e = pl.col(col).str.contains(needles[0], literal=True)
    for n in needles[1:]:
        e = e | pl.col(col).str.contains(n, literal=True)
    return e


def frequency() -> pl.Expr:
    """``method.frequency`` as a fraction in [0, 1]. Unparseable or absent is 0.0.

    Forms: ``n/m`` (and ``n//m``, which occurs), ``x%``, or a bare float.
    """
    f = pl.col("method.frequency").str.strip_chars()
    num = f.str.replace_all(r"/+", "/").str.split("/")
    parsed = (
        pl.when(f == "").then(0.0)
        .when(f.str.contains("/", literal=True))
        .then(num.list.get(0).cast(pl.Float64, strict=False)
              / num.list.get(1).cast(pl.Float64, strict=False))
        .when(f.str.ends_with("%") & (f.str.len_chars() > 1))
        .then(f.str.head(-1).cast(pl.Float64, strict=False) / 100.0)
        .otherwise(f.cast(pl.Float64, strict=False))
        .fill_nan(0.0).fill_null(0.0)
    )
    # A zero denominator gives inf, which clears every assay threshold below and scores the record
    # at the ceiling. `fill_nan` does not catch it. No chunk contains one today (measured
    # 2026-09-27: 0 of 192,753 rows, maximum frequency 1.0), but chunks are submitted, so the
    # unparseable-is-0.0 rule has to cover this form as well.
    return pl.when(parsed.is_finite()).then(parsed).otherwise(0.0)


def cell_count() -> pl.Expr:
    """The numerator of an ``n/m`` frequency. A percentage gives no cell count, so 0."""
    f = pl.col("method.frequency").str.strip_chars()
    return (
        pl.when(f.str.contains("/", literal=True))
        .then(f.str.replace_all(r"/+", "/").str.split("/").list.get(0)
              .cast(pl.Int64, strict=False))
        .otherwise(0).fill_null(0)
    )


def sequencing_score(freq: pl.Expr, count: pl.Expr) -> pl.Expr:
    """How much the sequence can be trusted: single cell > Sanger > amplicon depth."""
    single = _lower("method.singlecell")
    seq = _lower("method.sequencing")
    return (
        pl.when((single != "") & (single != "no")).then(3)
        .when(seq == "sanger").then(pl.when(count >= 2).then(3).otherwise(2))
        .when(seq == "amplicon-seq").then(pl.when(freq >= 0.01).then(3).otherwise(1))
        .otherwise(1)
    )


def _moderate_specificity(freq: pl.Expr) -> pl.Expr:
    """One point when the assay's own enrichment threshold is met. Culture is judged hardest."""
    m = "__ident"
    return (
        pl.when(_contains_any(m, _CULTURE_BASED)).then(pl.when(freq >= 0.5).then(1).otherwise(0))
        .when(_contains_any(m, _SORT_BASED)).then(pl.when(freq >= 0.05).then(1).otherwise(0))
        .when(_contains_any(m, _STIMULATION_BASED)).then(pl.when(freq >= 0.25).then(1).otherwise(0))
        .otherwise(0)
    )


def _high_specificity() -> pl.Expr:
    v = "__verif"
    return (
        pl.when(pl.col(v).str.contains("direct", literal=True)).then(3)
        .when(_contains_any(v, _STIMULATION_BASED)).then(2)
        .when(_contains_any(v, _SORT_BASED)).then(1)
        .otherwise(0)
    )


def row_score() -> pl.Expr:
    """The per-row score, before the per-signature maximum."""
    freq, count = pl.col("__freq"), pl.col("__count")
    seq = sequencing_score(freq, count)
    spec2 = _high_specificity()
    # A verified TCR was cloned, so its sequence is trusted regardless of how it was read.
    seq = pl.when(spec2 > 0).then(3).otherwise(seq)
    return (
        pl.when(pl.col("meta.structure.id") != "").then(3)
        .otherwise(pl.min_horizontal(seq, _moderate_specificity(freq) + spec2))
        .cast(pl.Int64)
    )


def add_score(df: pl.DataFrame) -> pl.DataFrame:
    """Add ``vdjdb.score``: the per-row score, maximised over :data:`SCORE_SIGNATURE`."""
    return (
        df.with_columns(
            frequency().alias("__freq"),
            cell_count().alias("__count"),
            _lower("method.identification").alias("__ident"),
            _lower("method.verification").alias("__verif"),
        )
        .with_columns(row_score().alias("__row_score"))
        .with_columns(pl.col("__row_score").max().over(SCORE_SIGNATURE).alias("vdjdb.score"))
        .drop("__freq", "__count", "__ident", "__verif", "__row_score")
    )
