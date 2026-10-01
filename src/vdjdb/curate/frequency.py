"""``method.frequency``: three independent columns, none derived from another (#696).

One free-text column holds at least three different things. Measured on the corpus: **63,532
non-blank cells over 7,111 distinct values**, of which 3,548 distinct are a ratio written `x/X`
(`12/163`, `52/27589`, `5719/33921`), others a bare float (`0.052115583`) and others a percentage
(`0.163%`, `18.06%`, `99%`).

**The count is the thing the submitter measured and the column VDJdb offers for it is a string.** A
clonotype supported by 3 reads is better evidence than one supported by 1, and that difference is
exactly what a frequency erases: at any realistic depth both are `~1/total` to several significant
figures, and `1/33921` against `3/33921` is 2.9e-5 against 8.8e-5, which a two-significant-figure
export makes the same number.

So ``method.frequency.count`` and ``method.frequency.total`` are parsed out where the submitted
string is an unambiguous ratio, and **the submitted string stays exactly as submitted**. Nothing is
derived in either direction:

* deriving the float from the pair would blank it on every study that reported only a float, which
  is where it is the only measurement there is;
* deriving the pair from the float is impossible.

Where all three are present they must **agree**, and that is a QC advisory rather than a rewrite -
a disagreement is a curation question about what the paper says, not something to round away.
"""
from __future__ import annotations

import polars as pl

#: An unambiguous count over a total. Anchored, so `1/2/3` and `12/163 of 200` are not parsed: a
#: partial read of a value nobody checked is worse than no read.
#:
#: ``//`` is accepted as well as ``/``. **297 records write the ratio with a doubled slash** - 106
#: distinct values, all of the form `1//13`, in `PMID_29150238.tsv` (239) and `PMID_11046006.tsv`
#: (58). A doubled separator carries no second
#: meaning, so reading it recovers 297 counts. The submitted string is **not** rewritten: only this
#: parse is tolerant, so `method.frequency` still shows what the curator typed and the repair, if
#: anyone wants one, stays a chunk edit with its own reason.
RATIO = r"^\s*(\d+)\s*//?\s*(\d+)\s*$"

#: A bare float, including the exponent form. **689 records of `PMID_28636589.tsv` over 17 distinct
#: values** are written `1e-04`, `2e-05` and so on; they are floats and the report must say so
#: rather than filing them as unreadable. They carry no count either way - a float never does.
FLOAT = r"^\s*\d*\.?\d+(?:[eE][-+]?\d+)?\s*$"

#: How far `count / total` may differ from a submitted float or percentage before the two are
#: reported as disagreeing. Relative, because the values span `1/33921` to `99%`; 1 % absorbs the
#: rounding in a two-significant-figure export without hiding a transcription error.
CONCORDANCE_TOLERANCE = 0.01


def split_frequency(df: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Add ``method.frequency.count`` and ``.total``; leave ``method.frequency`` untouched.

    Returns the frame and a report of what was parsed, by shape of the submitted value.
    """
    column = "method.frequency"
    empty = pl.DataFrame(schema={"shape": pl.String, "rows": pl.Int64})
    if column not in df.columns:
        return df, empty

    groups = pl.col(column).str.extract_groups(RATIO)
    # A submitted value wins over a parsed one: `method.frequency.count` and `.total` are chunk
    # columns, so a curator with a read count can write it rather than encode it in a string, and
    # this pass fills only what they left blank. Always trust the submission first.
    submitted = {c: (pl.col(c).cast(pl.Int64, strict=False) if c in df.columns else pl.lit(None))
                 for c in ("method.frequency.count", "method.frequency.total")}
    df = df.with_columns(
        submitted["method.frequency.count"]
            .fill_null(groups.struct[0].cast(pl.Int64)).alias("method.frequency.count"),
        submitted["method.frequency.total"]
            .fill_null(groups.struct[1].cast(pl.Int64)).alias("method.frequency.total"),
        submitted["method.frequency.count"].is_not_null().alias("__count.submitted"),
    )
    report = df.select(
        pl.when(pl.col("__count.submitted")).then(pl.lit("submitted count/total"))
          .when(pl.col(column) == "").then(pl.lit("blank"))
          .when(pl.col("method.frequency.count").is_not_null()).then(pl.lit("parsed count/total"))
          .when(pl.col(column).str.contains("%")).then(pl.lit("percentage"))
          .when(pl.col(column).str.contains(FLOAT)).then(pl.lit("float"))
          .otherwise(pl.lit("unparsed")).alias("shape"),
    ).group_by("shape").len().rename({"len": "rows"}).sort("rows", descending=True)
    return df.drop("__count.submitted"), report.select("shape", pl.col("rows").cast(pl.Int64))


def discordant() -> pl.Expr:
    """True where a submitted float or percentage disagrees with the count and total beside it.

    Live only because ``method.frequency.count`` and ``.total`` are chunk columns: a cell of
    ``method.frequency`` is one shape or the other, so before #696 no record could report all three
    and this could never fire. It reports and never rewrites - which of the three the paper supports
    is a curation question, not something to round away.
    """
    submitted = pl.col("method.frequency").str.strip_chars()
    as_fraction = (pl.when(submitted.str.ends_with("%"))
                     .then(submitted.str.strip_suffix("%").cast(pl.Float64, strict=False) / 100)
                     .otherwise(submitted.cast(pl.Float64, strict=False)))
    # Cast rather than assume: this expression is evaluated both on the typed record table and on
    # the all-string chunk frame `vdjdb qc` reads, where a division would raise.
    count = pl.col("method.frequency.count").cast(pl.Float64, strict=False)
    total = pl.col("method.frequency.total").cast(pl.Float64, strict=False)
    ratio = pl.when(total != 0).then(count / total).otherwise(None)
    return (as_fraction.is_not_null() & ratio.is_not_null()
            & ((as_fraction - ratio).abs() > ratio.abs() * CONCORDANCE_TOLERANCE))
