"""The dashboard's panels, in matplotlib.

Six of the dashboard's eight figures are ordinary geoms. They are drawn here instead of in R, in
the style the VDJdb papers already use: `~/vcs/manuscripts/2026-vdjdb-update` has no R at all --
every published figure is matplotlib with an Arial 7pt / 0.6pt-linewidth rcParams block and
`pdf.fonttype = 42` so the text stays editable. This module carries the same block, so a dashboard
panel and a paper panel are the same object.

**Two panels are deliberately not here.** The COVID and self-antigen alluvia need `ggalluvial` and
the TRBV-HLA chord needs `circlize`, and neither has a faithful Python counterpart -- plotly's
Sankey is a different object and every Python chord library draws a visibly different figure. Those
stay in R. Porting them "approximately" would change published figures to save a dependency, which
is the wrong trade.

Every function takes already-computed data and returns a `Figure`. The computation lives in
:func:`cumulative` and friends so it can be checked against the R that it replaces -- which it was:
the by-year numbers agree on all 408 (chain, metric, year) cells.

⚠ **Nothing here reads TCRvdb.**
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

#: Lifted from `2026-vdjdb-update/notebooks/vdjdb_timeline_figure.ipynb`. `fonttype 42` is what
#: keeps text as text in the PDF rather than outlines, which every journal asks for.
RC = {
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.dpi": 96,
    "savefig.bbox": "tight",
}

#: ColorBrewer Set1, the palette the R document uses for TCR chain. Written out rather than fetched
#: from a library so both renderers cannot drift (CLAUDE.md section 7 forbids hand-picked colours,
#: not named palettes).
SET1 = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
#: ColorBrewer PuBuGn, 4 classes -- the confidence-score fill.
PUBUGN4 = ["#f6eff7", "#bdc9e1", "#67a9cf", "#02818a"]


def style():
    """Apply :data:`RC`. Separate from the drawing so a caller can override before plotting."""
    import matplotlib

    matplotlib.rcParams.update(RC)


def cohort(legacy: Path) -> pl.DataFrame:
    """The slim table, which is what every panel below reads. One row per CDR3-antigen pair."""
    return pl.read_csv(legacy / "vdjdb.slim.txt", separator="\t", infer_schema_length=0,
                       quote_char=None)


def cumulative(legacy: Path, years: Path) -> pl.DataFrame:
    """``(chains, metric, year, total)`` -- distinct keys first seen at or before each year.

    The same quantity the R document plots, by the same method: a key's **first** year, tabulated
    and cumulated. The form it replaced cross-joined every (year, year) pair, roughly 7M rows.
    """
    full = pl.read_csv(legacy / "vdjdb_full.txt", separator="\t", infer_schema_length=0,
                       quote_char=None, truncate_ragged_lines=True)
    pub = (pl.read_csv(years, separator="\t")
             .select("reference.id", pl.col("year").cast(pl.Int32).alias("pub_date")))
    base = (full.filter(pl.col("species") != "MacacaMulatta")
            .select(
                "reference.id", "antigen.epitope",
                # `concat_str` with the nulls filled, NOT `+`: polars propagates a null through a
                # string concatenation, so one absent field made the whole key null and the rows
                # vanished at the `!= ""` filter. R's `paste` never does that, which is why the
                # port lost (TRA, tcr) and (TRB, tcr) entirely and got 25 other cells wrong.
                pl.concat_str(
                    [pl.col(c).fill_null("") for c in
                     ("v.alpha", "j.alpha", "cdr3.alpha", "v.beta", "j.beta", "cdr3.beta")],
                    separator=" ").alias("tcr"),
                pl.concat_str([pl.col("mhc.a").fill_null(""), pl.col("mhc.b").fill_null("")],
                              separator=" ").alias("mhc"),
                pl.when(pl.col("cdr3.alpha") != "")
                  .then(pl.when(pl.col("cdr3.beta") != "").then(pl.lit("paired"))
                        .otherwise(pl.lit("TRA")))
                  .otherwise(pl.lit("TRB")).alias("chains"))
            .unique()
            .join(pub, on="reference.id", how="inner"))
    out = []
    for metric, col in (("tcr", "tcr"), ("epi", "antigen.epitope"),
                        ("ref", "reference.id"), ("mhc", "mhc")):
        first = (base.filter(pl.col(col) != "")
                 .group_by("chains", pl.col(col).alias("key"))
                 .agg(pl.col("pub_date").min().alias("year")))
        out.append(first.group_by("chains", "year").len().rename({"len": "added"})
                        .with_columns(pl.lit(metric).alias("metric")))
    grid = pl.concat(out)
    span = pl.DataFrame({"year": sorted(pub["pub_date"].unique().to_list())})
    combos = (grid.select("chains").unique()
              .join(grid.select("metric").unique(), how="cross"))
    full_grid = (combos.join(span, how="cross")
                 .join(grid, on=["chains", "metric", "year"], how="left")
                 .with_columns(pl.col("added").fill_null(0))
                 .sort("chains", "metric", "year"))
    return full_grid.with_columns(
        pl.col("added").cum_sum().over("chains", "metric").alias("total"))


def by_year(cum: pl.DataFrame, annotations: pl.DataFrame | None = None):
    """The 2x2 cumulative grid: TCRs, epitopes, studies, MHC alleles."""
    import matplotlib.pyplot as plt

    style()
    titles = {"tcr": "Number of unique TCRs", "epi": "Number of unique epitopes",
              "ref": "Number of studies", "mhc": "Number of MHC alleles"}
    chains = ["TRA", "TRB", "paired"]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0))
    for ax, metric in zip(axes.ravel(), ("tcr", "epi", "ref", "mhc"), strict=True):
        d = cum.filter(pl.col("metric") == metric)
        for colour, chain in zip(SET1, chains, strict=False):
            c = d.filter(pl.col("chains") == chain).sort("year")
            ax.plot(c["year"], c["total"], "-o", color=colour, ms=2.5, lw=0.8, label=chain)
        ax.set_title(titles[metric], fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        # Every two years, as `scale_x_continuous(breaks = seq(1995, max, by = 2))` does. A denser
        # axis than matplotlib picks, and the one readers of this figure are used to.
        ax.set_xticks(range(1995, int(cum["year"].max()) + 1, 2))
        ax.tick_params(axis="x", rotation=90)
        if annotations is not None:
            _callouts(ax, d, annotations.filter(pl.col("panel") == metric))
    # One legend for the grid, at the bottom -- what `grid.arrange(..., mylegend)` produces in R.
    handles, labels = axes[0, 0].get_legend_handles_labels()
    # `rect` reserves the strip BEFORE tight_layout packs the axes; a bbox offset alone leaves the
    # legend sitting on top of the bottom row's tick labels.
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles, labels, title="TCR chain(s)", loc="lower center", ncols=len(labels),
               frameon=False, fontsize=6, title_fontsize=6)
    return fig


def _callouts(ax, panel_data: pl.DataFrame, marks: pl.DataFrame, float_frac: float = 0.05):
    """Segment to the series value in that year, label a fixed fraction of the panel above it.

    No coordinates in the table -- the same rule the R document uses, for the same reason (#460):
    the hardcoded positions it replaced were right for a 2022 database and wrong ever since.
    """
    if not marks.height or not panel_data.height:
        return
    top = panel_data["total"].max()
    for row in marks.iter_rows(named=True):
        at = panel_data.filter(pl.col("year") == row["year"])["total"]
        value = at.max() if at.len() else 0
        ax.plot([row["year"], row["year"]], [0, value], color="0.25", lw=0.3)
        ax.annotate(row["label"].replace("\\n", "\n"), (row["year"], value + float_frac * top),
                    ha="right", va="top", fontsize=6)
