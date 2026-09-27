"""The dashboard's panels, in matplotlib.

Five of the dashboard's eight figures are ordinary geoms. They are drawn here instead of in R, in
the style the VDJdb papers already use: `~/vcs/manuscripts/2026-vdjdb-update` has no R at all --
every published figure is matplotlib with an Arial 7pt / 0.6pt-linewidth rcParams block and
`pdf.fonttype = 42` so the text stays editable. This module repeats the same block, so a dashboard
panel and a paper panel are the same object.

Three panels are deliberately not here. The COVID and self-antigen alluvia need `ggalluvial`
and the TRBV-HLA chord needs `circlize`, and neither has a faithful Python counterpart -- plotly's
Sankey is a different object and every Python chord library draws a visibly different figure. Those
stay in R, because an approximate port would change published figures to save a dependency.

Every panel here was checked cell by cell against the R it replaces, on the current corpus:

=====================  ==============  =============
panel                  cells compared  differences
=====================  ==============  =============
by-year grid                      408              0
V-gene x MHC allele              2327              0
spectratype                        42              0
epitope length                     24              0
confidence score                   16              0
=====================  ==============  =============

Every function takes already-computed data and returns a `Figure`. The computation sits in
:func:`cumulative` and its siblings so it can be checked against the R it replaces: the by-year
numbers agree on all 408 (chain, metric, year) cells.

⚠ Nothing here reads TCRvdb.
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

    The same quantity the R document plots, by the same method: a key's first year, tabulated
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
                # string concatenation, so one absent field made the key null and the rows
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
        # axis than matplotlib picks, and the one readers of this figure expect.
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


def scores(cohort_df: pl.DataFrame):
    """Confidence-score distribution: grouped bars per (MHC class, chain), log records."""
    import matplotlib.pyplot as plt
    import numpy as np

    style()
    d = (cohort_df.filter(pl.col("species") == "HomoSapiens")
         .group_by("mhc.class", "gene", "vdjdb.score").len().rename({"len": "total"})
         .with_columns((pl.col("mhc.class") + " " + pl.col("gene")).alias("group"))
         .sort("group", "vdjdb.score"))
    groups = d["group"].unique(maintain_order=True).sort().to_list()
    levels = d["vdjdb.score"].unique().sort().to_list()
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    width = 0.8 / len(levels)
    for i, level in enumerate(levels):
        sub = d.filter(pl.col("vdjdb.score") == level)
        by_group = dict(zip(sub["group"], sub["total"], strict=True))
        x = np.arange(len(groups)) + (i - (len(levels) - 1) / 2) * width
        ax.bar(x, [by_group.get(g, 0) for g in groups], width=width,
               color=PUBUGN4[i % len(PUBUGN4)], edgecolor="black", linewidth=0.3, label=str(level))
    ax.set_yscale("log")
    ax.set_ylabel("Records")
    ax.set_xticks(np.arange(len(groups)), groups)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="VDJdb score", frameon=False, fontsize=6, title_fontsize=6,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncols=len(levels))
    fig.tight_layout()
    return fig


def nrd0(v) -> float:
    """R's default density bandwidth, ``bw.nrd0``: ``0.9 * min(sd, IQR/1.34) * n^(-1/5)``.

    Written out because the curve it produces is the one readers of this figure know, and because
    scipy has no equivalent -- its Scott and Silverman rules are both narrower, and on integer
    CDR3 lengths Scott's draws a comb of spikes that overshoots the bars.

    **The divisor is 1.34, not 1.349.** ``bw.nrd`` uses 1.349; ``bw.nrd0``, which is what
    ``geom_density`` defaults to, uses 1.34. Reading it off the documented formula instead of R's
    source put this 0.67 % out -- close enough to look correct. Checked against R's own output on
    two fixed vectors (``tests/unit/test_panels.py``).

    R falls back through sd, then |x[0]|, then 1 when the spread is zero; so does this.
    """
    import numpy as np

    n = len(v)
    sd = float(np.std(v, ddof=1))
    iqr = float(np.subtract(*np.percentile(v, [75, 25])))
    lo = min(sd, iqr / 1.34)
    if not lo:                      # R: `(lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)`
        lo = sd or abs(float(v[0])) or 1.0
    return 0.9 * lo * n ** (-0.2)


def spectratype(cohort_df: pl.DataFrame, *, lo: int = 5, hi: int = 25, adjust: float = 3.0):
    """CDR3 length distribution per chain, stacked by epitope and coloured by epitope length.

    The fill is a Spectral gradient over epitopes ordered by their own length, which is what
    `fct_reorder(epi_len) %>% as.integer()` does in the R -- the colour encodes epitope length,
    not identity, so the stack reads as "short epitopes at one end of the spectrum".

    The dotted overlay is a Gaussian KDE scaled to counts, using R's own bandwidth rule -- see
    :func:`nrd0`. scipy's default (Scott's rule) is far too narrow for integer lengths and draws a
    comb of spikes rather than a curve; measured on this data it overshot the bars by 2x.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colormaps
    from scipy.stats import gaussian_kde

    style()
    d = (cohort_df.filter(pl.col("species") == "HomoSapiens")
         .select("gene", "antigen.epitope",
                 pl.col("cdr3").str.len_chars().alias("len"),
                 pl.col("antigen.epitope").str.len_chars().alias("epi_len"))
         .filter(pl.col("len").is_between(lo, hi)))
    order = (d.select("antigen.epitope", "epi_len").unique()
             .sort("epi_len", "antigen.epitope")["antigen.epitope"].to_list())
    rank = {e: i for i, e in enumerate(order)}
    cmap = colormaps["Spectral"]
    genes = d["gene"].unique().sort().to_list()
    bins = np.arange(lo, hi + 2) - 0.5
    fig, axes = plt.subplots(1, len(genes), figsize=(6.0, 2.6), sharey=True)
    for ax, gene in zip(np.atleast_1d(axes), genes, strict=True):
        g = d.filter(pl.col("gene") == gene)
        bottom = np.zeros(len(bins) - 1)
        for epi, sub in g.group_by("antigen.epitope"):
            counts, _ = np.histogram(sub["len"].to_numpy(), bins=bins)
            ax.bar(bins[:-1] + 0.5, counts, bottom=bottom, width=1.0,
                   color=cmap(rank[epi[0]] / max(len(order) - 1, 1)), alpha=0.9, linewidth=0)
            bottom += counts
        xs = np.linspace(lo, hi, 200)
        v = g["len"].to_numpy().astype(float)
        # scipy's `bw_method` is a factor on the data's own standard deviation, so R's bandwidth
        # has to be divided by it to mean the same thing.
        kde = gaussian_kde(v, bw_method=(adjust * nrd0(v)) / v.std(ddof=1))
        ax.plot(xs, kde(xs) * g.height, ls=":", color="black", lw=0.6)
        ax.set_title(gene, fontsize=7)
        ax.set_xlabel("CDR3 length")
        ax.set_xticks(range(lo, hi + 1, 5))
        ax.spines[["top", "right"]].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("Records")
    fig.tight_layout()
    return fig


def v_hla(cohort_df: pl.DataFrame, *, min_records: int = 10) -> pl.DataFrame:
    """``(gene, mhc.class, mhc, v, records)`` for the V-gene x MHC-allele heatmap.

    Three columns are comma-separated lists and all three are exploded, so one slim row can land in
    several cells. ``records`` therefore counts distinct source rows, not exploded ones --
    `length(unique(id))` in the R, and the reason the row id is assigned before the explosion
    rather than after.

    Alleles are truncated at the first ``:`` (two-field resolution) and V genes at ``*``, then
    anything whose MHC or V total is below ``min_records`` is dropped, exactly as the R does.
    """
    d = (cohort_df.filter(pl.col("species") == "HomoSapiens")
         .with_row_index("id")
         .select("id", "gene", "mhc.class", "mhc.a", "mhc.b", "v.segm")
         .with_columns(pl.col("mhc.a", "mhc.b", "v.segm").str.split(","))
         # Pinned: Polars 2.0 flips the default. `str.split` on an empty cell yields `[""]`,
         # never an empty list, so this changes nothing here -- and pinning stops an upgrade from
         # moving a published figure with no error.
         .explode("mhc.a", empty_as_null=False)
         .explode("mhc.b", empty_as_null=False)
         .explode("v.segm", empty_as_null=False)
         .with_columns(
             pl.col("mhc.a").str.split(":").list.first().alias("a"),
             pl.col("mhc.b").str.split(":").list.first().alias("b"),
             pl.col("v.segm").str.split("*").list.first().alias("v"))
         .with_columns((pl.col("a") + " / " + pl.col("b")).str.replace_all("HLA-", "",
                                                                          literal=True)
                       .alias("mhc")))
    cells = (d.group_by("gene", "mhc.class", "mhc", "v")
             .agg(pl.col("id").n_unique().alias("records")))
    by_mhc = cells.group_by("mhc.class", "mhc").agg(pl.col("records").sum().alias("mhc_total"))
    by_v = cells.group_by("gene", "v").agg(pl.col("records").sum().alias("v_total"))
    return (cells.join(by_mhc, on=["mhc.class", "mhc"])
            .join(by_v, on=["gene", "v"])
            .filter((pl.col("mhc_total") >= min_records) & (pl.col("v_total") >= min_records)))


def v_hla_heatmap(cells: pl.DataFrame, *, cap: int = 1000):
    """V gene x MHC allele, tiles on a log colour scale, faceted gene x MHC class.

    ``facet_grid(scales = "free", space = "free")`` has no matplotlib equivalent, so the panels are
    laid out on a GridSpec whose row and column ratios are the category counts -- which is what
    "free space" means: a tile is the same size in every panel.

    Axis order is ``fct_reorder(records)``, i.e. by the **median** records of each level, not the
    sum. Counts are capped at ``cap`` before colouring (`pmin(records, 1000)` in the R), so a
    handful of very large cells cannot flatten the rest of the scale.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colormaps
    from matplotlib.colors import LogNorm

    style()
    genes = cells["gene"].unique().sort().to_list()
    classes = cells["mhc.class"].unique().sort().to_list()
    # One shared axis order per column / per row, or tiles would not line up across facets.
    x_order = {c: (cells.filter(pl.col("mhc.class") == c).group_by("mhc")
                   .agg(pl.col("records").median().alias("m")).sort("m")["mhc"].to_list())
               for c in classes}
    y_order = {g: (cells.filter(pl.col("gene") == g).group_by("v")
                   .agg(pl.col("records").median().alias("m")).sort("m")["v"].to_list())
               for g in genes}

    fig = plt.figure(figsize=(6.0, 10.0))
    gs = fig.add_gridspec(len(genes), len(classes),
                          width_ratios=[max(len(x_order[c]), 1) for c in classes],
                          height_ratios=[max(len(y_order[g]), 1) for g in genes],
                          hspace=0.08, wspace=0.06)
    norm = LogNorm(vmin=1, vmax=cap)
    cmap = colormaps["PuBuGn"]
    mesh = None
    for r, gene in enumerate(genes):
        for c, klass in enumerate(classes):
            ax = fig.add_subplot(gs[r, c])
            xs, ys = x_order[klass], y_order[gene]
            grid = np.full((len(ys), len(xs)), np.nan)
            xi = {v: i for i, v in enumerate(xs)}
            yi = {v: i for i, v in enumerate(ys)}
            for row in cells.filter((pl.col("gene") == gene)
                                    & (pl.col("mhc.class") == klass)).iter_rows(named=True):
                grid[yi[row["v"]], xi[row["mhc"]]] = min(row["records"], cap)
            mesh = ax.pcolormesh(np.arange(len(xs) + 1), np.arange(len(ys) + 1), grid,
                                 cmap=cmap, norm=norm, edgecolors="none")
            # Tick labels on the OUTSIDE edges only. Labelling every panel puts the right
            # column's y-axis text on top of the left column's tiles, and repeats the allele names
            # on both rows, which `facet_grid` avoids by construction.
            if r == len(genes) - 1:
                ax.set_xticks(np.arange(len(xs)) + 0.5, xs, rotation=90, fontsize=5)
            else:
                ax.set_xticks([])
            if c == 0:
                ax.set_yticks(np.arange(len(ys)) + 0.5, ys, fontsize=5)
            else:
                ax.set_yticks([])
            if r == 0:
                ax.set_title(klass, fontsize=7)
            if c == len(classes) - 1:
                # The right-hand strip `facet_grid(gene ~ .)` draws.
                ax.text(1.01, 0.5, gene, transform=ax.transAxes, rotation=270,
                        va="center", ha="left", fontsize=7)
            ax.spines[["top", "right"]].set_visible(False)
    if mesh is not None:
        fig.colorbar(mesh, ax=fig.axes, label="Records", fraction=0.025, pad=0.02,
                     ticks=[1, 10, 100, 1000])
    return fig


def epitope_length(cohort_df: pl.DataFrame):
    """Epitope length per MHC class, stacked by CDR3 and coloured by CDR3 length.

    The R stacks one bar segment per distinct CDR3 -- upwards of a hundred thousand of them --
    which is why the panel reads as fine horizontal striations rather than flat blocks. Drawing
    that many matplotlib artists is not viable, so each bin is drawn as a one-pixel-wide **image**
    whose rows are the records in the same rank order. Same picture, O(bins) draw calls instead of
    O(records).

    Rank is ``fct_reorder(cdr3, cdr3_len)``: CDR3s ordered by their own length, ties alphabetical.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colormaps

    style()
    d = (cohort_df.filter(pl.col("species") == "HomoSapiens")
         .select("mhc.class", "cdr3",
                 pl.col("antigen.epitope").str.len_chars().alias("epi_len"),
                 pl.col("cdr3").str.len_chars().alias("cdr3_len"))
         .filter(pl.col("epi_len") > 0))
    order = (d.select("cdr3", "cdr3_len").unique()
             .sort("cdr3_len", "cdr3")["cdr3"].to_list())
    rank = {c: i for i, c in enumerate(order)}
    cmap = colormaps["Spectral"]
    classes = d["mhc.class"].unique().sort().to_list()
    fig, axes = plt.subplots(1, len(classes), figsize=(6.5, 2.8))
    for ax, klass in zip(np.atleast_1d(axes), classes, strict=True):
        g = d.filter(pl.col("mhc.class") == klass)
        lengths = sorted(g["epi_len"].unique().to_list())
        for n in lengths:
            ranks = np.sort(np.array([rank[c] for c in
                                      g.filter(pl.col("epi_len") == n)["cdr3"].to_list()]))
            if not ranks.size:
                continue
            col = cmap(ranks / max(len(order) - 1, 1))[:, None, :]
            ax.imshow(col, origin="lower", aspect="auto", interpolation="nearest",
                      extent=(n - 0.5, n + 0.5, 0, ranks.size))
        ax.set_xlim(min(lengths) - 0.6, max(lengths) + 0.6)
        ax.set_ylim(0, None)
        ax.autoscale(axis="y")
        ax.set_title(klass, fontsize=7)
        ax.set_xlabel("Epitope length")
        ax.set_xticks(lengths)
        ax.tick_params(axis="x", labelsize=5)
        ax.spines[["top", "right"]].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("Records")
    fig.tight_layout()
    return fig
