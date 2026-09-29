"""The interactive dashboard: that it is built from the verified data and that it is one file."""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

pytest.importorskip("plotly", reason="needs the `summary` extra")
sys.path.insert(0, "summary")

from vdjdb.summary import interactive


def _legacy(tmp: Path) -> Path:
    """A two-record legacy projection, the same fixture shape `test_panels.py` builds."""
    d = tmp / "legacy"
    d.mkdir()
    slim = pl.DataFrame({
        "gene": ["TRA", "TRB"], "cdr3": ["CAVRDSNYQLIW", "CASSIRSSYEQYF"],
        "species": ["HomoSapiens", "HomoSapiens"],
        "antigen.epitope": ["GILGFVFTL", "NLVPMVATV"],
        "antigen.gene": ["M", "pp65"], "antigen.species": ["InfluenzaA", "CMV"],
        "complex.id": ["0", "0"], "v.segm": ["TRAV12-2*01", "TRBV19*01"],
        "j.segm": ["TRAJ33*01", "TRBJ2-7*01"], "mhc.a": ["HLA-A*02:01", "HLA-A*02:01"],
        "mhc.b": ["B2M", "B2M"], "mhc.class": ["MHCI", "MHCI"],
        "reference.id": ["PMID:1", "PMID:2"], "vdjdb.score": ["1", "2"],
        "TCR_hash": ["a", "b"], "j.start": ["7", "8"], "v.end": ["3", "4"]})
    slim.write_csv(d / "vdjdb.slim.txt", separator="\t")
    full = pl.DataFrame({
        "reference.id": ["PMID:1", "PMID:2"], "antigen.epitope": ["GILGFVFTL", "NLVPMVATV"],
        "species": ["HomoSapiens", "HomoSapiens"],
        "v.alpha": ["TRAV12-2*01", ""], "j.alpha": ["TRAJ33*01", ""],
        "cdr3.alpha": ["CAVRDSNYQLIW", ""], "v.beta": ["", "TRBV19*01"],
        "j.beta": ["", "TRBJ2-7*01"], "cdr3.beta": ["", "CASSIRSSYEQYF"],
        "mhc.a": ["HLA-A*02:01", "HLA-A*02:01"], "mhc.b": ["B2M", "B2M"]})
    full.write_csv(d / "vdjdb_full.txt", separator="\t")
    return d


@pytest.fixture
def years(tmp_path: Path) -> Path:
    p = tmp_path / "reference_years.tsv"
    pl.DataFrame({"reference.id": ["PMID:1", "PMID:2"], "year": [2015, 2020]}).write_csv(
        p, separator="\t")
    return p


def test_the_palettes_match_the_static_dashboard_exactly():
    """Both dashboards declare the same ColorBrewer literals, so a chain is one colour everywhere.

    They are declared twice on purpose - importing `panels` pulls matplotlib, which is in an extra -
    and this is what stops the copies drifting.
    """
    import panels

    assert interactive.SET1 == panels.SET1
    assert interactive.PUBUGN4 == panels.PUBUGN4


def test_it_writes_one_self_contained_file_that_loads_plotly_from_a_cdn(tmp_path, years):
    """One `<script src>` for the whole page, not one per panel, and no bundled plotly.js.

    Bundling would make the file ~4 MB; six panels each requesting the CDN would be five wasted
    round trips. GitHub Pages serves the result as an ordinary static asset.
    """
    out = tmp_path / "interactive.html"
    interactive.build(_legacy(tmp_path), years, out)
    html = out.read_text()
    assert html.count("<script") - html.count("cdn.plot.ly") >= 1
    assert html.count("cdn.plot.ly") == 1, "plotly.js is requested once for the page"
    assert "plotly.min.js" not in html.split("cdn.plot.ly")[0], "plotly.js must not be bundled"
    assert len(html) < 2_000_000


def test_every_declared_section_reaches_the_page(tmp_path, years):
    """`SECTIONS` is the page order; a panel declared and not rendered would be a silent hole."""
    out = tmp_path / "interactive.html"
    html = interactive.build(_legacy(tmp_path), years, out).read_text()
    for key, heading, _ in interactive.SECTIONS:
        assert f'id="{key}"' in html, key
        assert heading in html, heading
    assert html.count("plotly-graph-div") == len(interactive.SECTIONS)


def test_the_header_counts_come_from_the_cohort(tmp_path, years):
    out = tmp_path / "interactive.html"
    html = interactive.build(_legacy(tmp_path), years, out).read_text()
    assert "2 records over 2 publications and 2 epitopes" in html


def test_the_heatmap_hover_reports_the_uncapped_count(tmp_path, years):
    """Capping is for the colour scale only. A reader hovering a capped cell must see the real
    number, or the cap becomes something they have to know about to read the panel."""
    cells = pl.DataFrame({"gene": ["TRB"], "mhc.class": ["MHCI"], "mhc": ["A*02"],
                          "v": ["TRBV19"], "records": [5000], "mhc_total": [5000],
                          "v_total": [5000]})
    fig = interactive.v_hla(cells, cap=1000)
    trace = fig.data[0]
    assert trace.z[0][0] == 1000, "the colour scale sees the capped value"
    assert trace.customdata[0][0] == 5000, "the hover sees the real value"


def test_a_missing_cell_is_a_gap_rather_than_a_zero(tmp_path, years):
    """A V/MHC pair with no records did not occur; drawing it as 0 claims it was looked for."""
    cells = pl.DataFrame({"gene": ["TRB", "TRB"], "mhc.class": ["MHCI", "MHCI"],
                          "mhc": ["A*02", "B*07"], "v": ["TRBV19", "TRBV20"],
                          "records": [10, 20], "mhc_total": [10, 20], "v_total": [10, 20]})
    fig = interactive.v_hla(cells)
    assert fig.data[0].hoverongaps is False
    assert any(v is None for line in fig.data[0].z for v in line)


def test_a_chain_keeps_one_colour_and_one_legend_entry_across_the_metrics(tmp_path, years):
    """`legendgroup` is what makes one click isolate a chain in all four metrics at once."""
    cum = pl.DataFrame({
        "chains": ["TRA"] * 4 + ["TRB"] * 4,
        "metric": ["tcr", "tcr", "epi", "epi"] * 2,
        "year": [2015, 2020] * 4, "added": [1, 1] * 4, "total": [1, 2] * 4})
    fig = interactive.by_year(cum)
    shown = [t for t in fig.data if t.showlegend]
    assert len(shown) == 2, "one legend entry per chain, not one per chain and metric"
    for trace in fig.data:
        assert trace.legendgroup == trace.name
        assert trace.line.color == interactive.CHAIN_COLOUR[trace.name]
    assert fig.layout.hovermode == "x unified"
