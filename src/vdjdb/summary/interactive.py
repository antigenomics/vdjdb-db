"""The interactive dashboard: one self-contained HTML file, served from GitHub Pages.

An **additive** artifact. The R document stays the shipped dashboard - `vdjdb-web` injects its
fragment into `/overview` and `summary/check_summary.py` gates it - and nothing here feeds that path.
This is the version you can interrogate: toggle a series, read an exact count off a cell, sort a
table. The static panels cannot do any of that, and three of them are hard to read because of it.

Where interaction actually buys something, which is the only reason to build a second dashboard:

- **V-gene x MHC allele.** 8 pt labels on a 2,327-cell grid are unreadable in print; hover gives the
  gene, the allele and the exact record count.
- **Cumulative by year.** Four metrics x three chain configurations on one set of axes. Clicking a
  chain in the legend isolates it; a unified hover reads all three at one year.
- **The HLA table.** Long, and worth sorting by a column other than the one it was written in.

The data comes from :mod:`summary.panels`, not from a second set of queries. Those functions were
checked cell by cell against the R they replace - 408 by-year cells, 2,327 V x MHC cells, 42
spectratype, 24 epitope length, 16 score, zero differences - so reusing them is what keeps this
dashboard and the static one from disagreeing about a number. Adding a query here would throw that
away.

`plotly.js` is loaded from a CDN (`include_plotlyjs="cdn"`), so the file is ~200 KB rather than ~4 MB
and GitHub Pages serves it as an ordinary static asset. No server, no Quarto, no Jupyter, no build
step beyond running this.

⚠ Nothing here reads TCRvdb, and nothing here reads a background: the dashboard is derived statistics
over `chunks/` only.
"""
from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType

import polars as pl

from .render import SUMMARY


def _panels() -> ModuleType:
    """``summary/panels.py``, which is a working script beside the Rmd rather than shipped code.

    `pyproject.toml` packages `src/vdjdb` only, so `panels` is not importable as
    `vdjdb.summary.panels`. Loaded the way ``tests/unit/test_panels.py`` loads it, and lazily,
    because it imports matplotlib and that lives in the `summary` extra.
    """
    if str(SUMMARY) not in sys.path:
        sys.path.insert(0, str(SUMMARY))
    import panels

    return panels

#: `theme_classic()` plus `axis.line = element_line(linewidth = 0.3)`, which is what the R document
#: sets: white panel, no gridlines, a black axis line on the left and bottom only, ticks outside.
#: Declared once and applied to every figure, so the two dashboards look like one product.
TEMPLATE = "vdjdb_classic"

#: Metric key -> the label the R document uses. Same order as its facet strip.
METRICS: dict[str, str] = {"tcr": "TCR chains", "epi": "Epitopes",
                           "mhc": "MHC alleles", "ref": "Publications"}

#: ColorBrewer Set1 and PuBuGn, the two palettes the R document and ``summary/panels.py`` use.
#: Repeated here rather than imported from `panels`, which pulls matplotlib at import time; the
#: duplication is asserted away in ``tests/unit/test_interactive.py``.
SET1 = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
PUBUGN4 = ["#f6eff7", "#bdc9e1", "#67a9cf", "#02818a"]

#: Chain configuration -> colour, so a chain is the same colour in both dashboards.
CHAIN_COLOUR: dict[str, str] = dict(zip(("paired", "TRA", "TRB"), SET1, strict=False))


def _register_template() -> None:
    """Install :data:`TEMPLATE` into plotly's registry."""
    import plotly.graph_objects as go
    import plotly.io as pio

    axis = {"showline": True, "linecolor": "black", "linewidth": 0.8, "mirror": False,
            "showgrid": False, "zeroline": False, "ticks": "outside", "tickcolor": "black",
            "automargin": True}
    pio.templates[TEMPLATE] = go.layout.Template(layout={
        "font": {"family": "Arial, Helvetica, sans-serif", "size": 12, "color": "#1b1b1b"},
        "paper_bgcolor": "white", "plot_bgcolor": "white",
        "xaxis": axis, "yaxis": axis,
        "colorway": SET1,
        "margin": {"l": 60, "r": 20, "t": 40, "b": 50},
        "hoverlabel": {"font": {"family": "Arial, Helvetica, sans-serif", "size": 12}},
        "legend": {"bgcolor": "rgba(0,0,0,0)", "borderwidth": 0},
    })
    pio.templates.default = TEMPLATE


def _nothing_to_show(message: str):
    """A valid, empty figure carrying an explanation.

    Every panel is built from a filtered subset, so any of them can legitimately be empty - and a
    figure that raises takes the whole page with it rather than losing one panel.
    """
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, xref="paper", yref="paper", x=0.5, y=0.5)
    fig.update_layout(height=140, xaxis={"visible": False}, yaxis={"visible": False})
    return fig


def by_year(cum: pl.DataFrame):
    """Cumulative distinct keys by year: one subplot per metric, one trace per chain configuration.

    The interaction is the reason this panel is here. `legendgroup` ties a chain's four traces
    together so one click isolates it across every metric at once, and `hovermode="x unified"` reads
    all three chain configurations at a single year rather than whichever line the cursor is nearest.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    keys = [m for m in METRICS if m in set(cum["metric"])]
    if not keys:
        return _nothing_to_show("No publication years resolved, so nothing can be plotted by year.")
    fig = make_subplots(rows=1, cols=len(keys), shared_xaxes=True,
                        subplot_titles=[METRICS[m] for m in keys], horizontal_spacing=0.06)
    for col, metric in enumerate(keys, start=1):
        for chains in sorted(cum["chains"].unique().to_list()):
            d = (cum.filter((pl.col("metric") == metric) & (pl.col("chains") == chains))
                 .sort("year"))
            if d.is_empty():
                continue
            fig.add_trace(go.Scatter(
                x=d["year"].to_list(), y=d["total"].to_list(), name=chains,
                legendgroup=chains, showlegend=col == 1, mode="lines",
                line={"width": 2, "color": CHAIN_COLOUR.get(chains)},
                hovertemplate=f"{chains}: %{{y:,}}<extra></extra>"), row=1, col=col)
        fig.update_xaxes(title_text="Year", row=1, col=col)
    fig.update_yaxes(title_text="Cumulative records", row=1, col=1)
    # No figure title: the section heading above the panel already carries it, and a centred title
    # over a subplot grid lands on top of the facet strip.
    fig.update_layout(height=380, hovermode="x unified",
                      legend={"orientation": "h", "y": -0.22, "x": 0.5, "xanchor": "center"})
    return fig


def v_hla(cells: pl.DataFrame, *, cap: int = 1000):
    """V gene x MHC allele, as one readable heatmap per (chain, MHC class).

    The static panel labels 2,343 cells at 8 pt, which is the one place the print dashboard is simply
    unreadable. Counts are capped for the colour scale exactly as the static panel caps them, because
    a handful of very large cells otherwise flatten everything else - but the hover reports the
    **uncapped** count, so capping stops being something the reader has to know about.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    facets = (cells.select("gene", "mhc.class").unique()
              .sort("gene", "mhc.class").rows())
    if not facets:
        # No cell clears the threshold. `make_subplots(rows=0)` raises, so a corpus small enough to
        # have no qualifying pair would take the whole page down rather than showing an empty panel.
        return _nothing_to_show("No V gene and MHC allele pair has 10 or more records.")
    fig = make_subplots(rows=len(facets), cols=1, vertical_spacing=0.06,
                        subplot_titles=[f"{g} x {c}" for g, c in facets])
    for row, (gene, mhc_class) in enumerate(facets, start=1):
        d = cells.filter((pl.col("gene") == gene) & (pl.col("mhc.class") == mhc_class))
        alleles = sorted(d["mhc"].unique().to_list())
        genes = sorted(d["v"].unique().to_list())
        counts = {(r["v"], r["mhc"]): r["records"] for r in d.iter_rows(named=True)}
        raw = [[counts.get((v, a)) for a in alleles] for v in genes]
        fig.add_trace(go.Heatmap(
            z=[[None if n is None else min(n, cap) for n in line] for line in raw],
            x=alleles, y=genes, customdata=raw, colorscale="PuBuGn", coloraxis="coloraxis",
            hovertemplate="%{y} x %{x}<br>%{customdata:,} records<extra></extra>",
            hoverongaps=False), row=row, col=1)
        fig.update_xaxes(tickangle=-90, tickfont={"size": 9}, row=row, col=1)
        fig.update_yaxes(tickfont={"size": 9}, row=row, col=1)
    heights = [max(160, 13 * cells.filter((pl.col("gene") == g)
                                          & (pl.col("mhc.class") == c))["v"].n_unique())
               for g, c in facets]
    fig.update_layout(
        height=sum(heights) + 120 * len(facets),
        coloraxis={"colorscale": "PuBuGn",
                   "colorbar": {"title": f"Records<br>(capped {cap:,})", "thickness": 12}})
    return fig


def spectratype(cohort: pl.DataFrame, *, lo: int = 5, hi: int = 25):
    """CDR3 length per chain. Counts, not the stacked epitope gradient the static panel draws.

    Stated rather than smuggled: the static panel stacks one stratum per epitope so the fill encodes
    epitope length. That gradient is a print device - here the same information is a hover and a
    legend, and drawing 1,946 strata in the browser would make the page slow for nothing.
    """
    import plotly.graph_objects as go

    d = (cohort.filter(pl.col("species") == "HomoSapiens")
         .select("gene", pl.col("cdr3").str.len_chars().alias("len"))
         .filter(pl.col("len").is_between(lo, hi))
         .group_by("gene", "len").agg(pl.len().alias("records"))
         .sort("gene", "len"))
    fig = go.Figure()
    for i, gene in enumerate(sorted(d["gene"].unique().to_list())):
        sub = d.filter(pl.col("gene") == gene).sort("len")
        fig.add_trace(go.Bar(x=sub["len"].to_list(), y=sub["records"].to_list(), name=gene,
                             marker={"color": SET1[i % len(SET1)]},
                             hovertemplate=f"{gene}, CDR3 length %{{x}}: %{{y:,}}<extra></extra>"))
    fig.update_layout(height=360, barmode="group",
                      xaxis={"title": "CDR3 length", "dtick": 5},
                      yaxis={"title": "Records"})
    return fig


def epitope_length(cohort: pl.DataFrame):
    """Epitope length per MHC class, with the distinct-CDR3 count on the hover.

    The static panel encodes distinct CDR3s as striations, one per CDR3 - which is what makes that
    panel 137.82 s of a 152 s R render (#648). Here it is a number in the tooltip.
    """
    import plotly.graph_objects as go

    d = (cohort.filter(pl.col("species") == "HomoSapiens")
         .select("mhc.class", "cdr3", pl.col("antigen.epitope").str.len_chars().alias("epi_len"))
         .filter(pl.col("epi_len") > 0)
         .group_by("mhc.class", "epi_len")
         .agg(pl.len().alias("records"), pl.col("cdr3").n_unique().alias("cdr3s"))
         .sort("mhc.class", "epi_len"))
    fig = go.Figure()
    for i, cls in enumerate(sorted(d["mhc.class"].unique().to_list())):
        sub = d.filter(pl.col("mhc.class") == cls)
        fig.add_trace(go.Bar(
            x=sub["epi_len"].to_list(), y=sub["records"].to_list(), name=cls,
            marker={"color": SET1[i % len(SET1)]},
            customdata=sub["cdr3s"].to_list(),
            hovertemplate=(f"{cls}, epitope length %{{x}}<br>%{{y:,}} records"
                           "<br>%{customdata:,} distinct CDR3<extra></extra>")))
    fig.update_layout(height=360, barmode="group",
                      xaxis={"title": "Epitope length", "dtick": 1},
                      yaxis={"title": "Records"})
    return fig


def scores(cohort: pl.DataFrame):
    """Confidence score by MHC class and chain, log records - the same grouping as the static panel."""
    import plotly.graph_objects as go

    d = (cohort.filter(pl.col("species") == "HomoSapiens")
         .group_by("mhc.class", "gene", "vdjdb.score").agg(pl.len().alias("records"))
         .with_columns((pl.col("mhc.class") + " " + pl.col("gene")).alias("group"))
         .sort("group", "vdjdb.score"))
    groups = sorted(d["group"].unique().to_list())
    fig = go.Figure()
    for i, level in enumerate(sorted(d["vdjdb.score"].unique().to_list())):
        sub = d.filter(pl.col("vdjdb.score") == level)
        got = dict(zip(sub["group"].to_list(), sub["records"].to_list(), strict=True))
        fig.add_trace(go.Bar(
            x=groups, y=[got.get(g, 0) for g in groups], name=str(level),
            marker={"color": PUBUGN4[i % len(PUBUGN4)],
                    "line": {"color": "black", "width": 0.5}},
            hovertemplate=f"score {level}, %{{x}}: %{{y:,}}<extra></extra>"))
    fig.update_layout(height=360, barmode="group",
                      yaxis={"title": "Records", "type": "log"},
                      legend={"title": "VDJdb score"})
    return fig


def hla_table(cohort: pl.DataFrame, *, top: int = 40):
    """The MHC allele table, sortable in the browser. `go.Table` and not HTML, so it needs no CSS."""
    import plotly.graph_objects as go

    d = (cohort.filter(pl.col("mhc.a") != "")
         .group_by("mhc.a")
         .agg(pl.len().alias("records"),
              pl.col("antigen.epitope").n_unique().alias("epitopes"),
              pl.col("cdr3").n_unique().alias("cdr3s"))
         .sort("records", descending=True).head(top))
    fig = go.Figure(go.Table(
        header={"values": ["<b>MHC allele</b>", "<b>Records</b>", "<b>Epitopes</b>",
                           "<b>Distinct CDR3</b>"],
                "fill_color": PUBUGN4[1], "align": "left"},
        cells={"values": [d["mhc.a"].to_list(),
                          [f"{v:,}" for v in d["records"]],
                          [f"{v:,}" for v in d["epitopes"]],
                          [f"{v:,}" for v in d["cdr3s"]]],
               "align": "left", "height": 22}))
    fig.update_layout(height=min(900, 40 + 22 * (len(d) + 1)), margin={"t": 10})
    return fig


#: Panel id -> (heading, one-line note). The order is the page order.
SECTIONS: tuple[tuple[str, str, str], ...] = (
    ("by-year", "Growth by year",
     "Distinct TCR chains, epitopes, MHC alleles and publications first reported at or before each "
     "year. Click a chain configuration in the legend to isolate it across all four metrics."),
    ("v-hla", "V gene and MHC allele",
     "One heatmap per chain and MHC class, cells with at least 10 records. The colour scale is "
     "capped so a few very large cells do not flatten the rest; the hover reports the uncapped "
     "count."),
    ("spectratype", "CDR3 length", "Homo sapiens, both chains."),
    ("epitope-length", "Epitope length",
     "By MHC class. The hover carries the number of distinct CDR3 behind each bar."),
    ("scores", "Confidence score", "Log record count, by MHC class and chain."),
    ("hla-table", "MHC alleles", "Click a column header to sort."),
)

_PAGE = """<!DOCTYPE html>
<meta charset="utf-8">
<title>VDJdb summary - interactive</title>
<style>
  body {{ font: 15px/1.55 Arial, Helvetica, sans-serif; color: #1b1b1b; margin: 0 auto;
         max-width: 1180px; padding: 1.5rem 1.25rem 4rem; }}
  h1 {{ font-size: 1.7rem; margin: 0 0 .25rem; }}
  h2 {{ font-size: 1.15rem; margin: 2.25rem 0 .2rem; }}
  p.note {{ color: #555; margin: .2rem 0 .6rem; max-width: 70ch; }}
  p.meta {{ color: #555; font-size: .9rem; }}
  .panel {{ border: 1px solid #e3e3e3; border-radius: 4px; padding: .25rem; }}
  a {{ color: #253494; }}
</style>
<h1>VDJdb summary</h1>
<p class="meta">{records} records over {publications} publications and {epitopes} epitopes.
Built from <code>chunks/</code> by <code>vdjdb summary --interactive</code>. The printed dashboard is
the <a href="dashboard.html">static version</a>; this one is for reading exact values.</p>
{body}
<p class="meta">Derived statistics only. No background and no third-party database is included.</p>
"""


def build(legacy: Path, years: Path, out: Path, *, cap: int = 1000) -> Path:
    """Write the one-file interactive dashboard and return its path.

    Reads only the legacy projection of this build, which is what the static dashboard reads too, so
    the two cannot be looking at different data.
    """
    from plotly.io import to_html

    _register_template()
    panels = _panels()
    cohort = panels.cohort(legacy)
    cum = panels.cumulative(legacy, years)
    figures = {
        "by-year": by_year(cum),
        "v-hla": v_hla(panels.v_hla(cohort), cap=cap),
        "spectratype": spectratype(cohort),
        "epitope-length": epitope_length(cohort),
        "scores": scores(cohort),
        "hla-table": hla_table(cohort),
    }
    blocks = []
    for key, heading, note in SECTIONS:
        # `include_plotlyjs` on the first figure only, or the page carries one <script src> per
        # panel: harmless but six identical requests.
        html = to_html(figures[key], include_plotlyjs="cdn" if not blocks else False,
                       full_html=False, config={"displaylogo": False, "responsive": True})
        blocks.append(f'<h2 id="{key}">{heading}</h2>\n<p class="note">{note}</p>\n'
                      f'<div class="panel">{html}</div>')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(_PAGE.format(
        records=f"{cohort.height:,}",
        publications=f"{cohort['reference.id'].n_unique():,}",
        epitopes=f"{cohort['antigen.epitope'].n_unique():,}",
        body="\n".join(blocks)))
    return out
