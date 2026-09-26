"""The matplotlib panels: the computation they draw, and the null trap that broke it once.

The by-year computation was verified against the R it replaces on the real corpus -- 408
(chain, metric, year) cells, zero differences. That check needs R and the full build, so it is not
in the suite. What is here is the property that made the first port silently wrong, plus enough of
a smoke test that the figure cannot stop being drawn without something failing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

sys.path.insert(0, "summary")

matplotlib = pytest.importorskip("matplotlib", reason="needs the `summary` extra")
matplotlib.use("Agg")
import panels  # noqa: E402  - must follow the backend selection


def _legacy(tmp_path: Path) -> Path:
    """A two-record legacy projection: one paired, one beta-only with EMPTY alpha fields."""
    d = tmp_path / "legacy"
    d.mkdir()
    cols = ["species", "reference.id", "antigen.epitope", "v.alpha", "j.alpha", "cdr3.alpha",
            "v.beta", "j.beta", "cdr3.beta", "mhc.a", "mhc.b"]
    rows = [
        ["HomoSapiens", "PMID:1", "GILGFVFTL", "TRAV1", "TRAJ1", "CAAA",
         "TRBV1", "TRBJ1", "CASSA", "HLA-A*02:01", "B2M"],
        # The beta-only record: every alpha column empty. This is the row that vanished.
        ["HomoSapiens", "PMID:2", "NLVPMVATV", "", "", "",
         "TRBV2", "TRBJ2", "CASSB", "HLA-A*02:01", "B2M"],
    ]
    (d / "vdjdb_full.txt").write_text(
        "\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return d


@pytest.fixture
def years(tmp_path):
    p = tmp_path / "years.tsv"
    p.write_text("reference.id\tyear\tsource\nPMID:1\t2020\tpubmed\nPMID:2\t2021\tpubmed\n")
    return p


def test_a_chain_with_empty_columns_still_gets_a_key(tmp_path, years):
    """polars propagates nulls through `+`; R's `paste` does not.

    Building the TCR key with `+` made it null for any record with an absent field, the `!= ""`
    filter dropped those rows, and two whole (chain, metric) series disappeared from the figure
    with 25 other cells wrong. Measured against the R before the fix.
    """
    cum = panels.cumulative(_legacy(tmp_path), years)
    got = {(r["chains"], r["metric"]) for r in cum.iter_rows(named=True)}
    assert ("TRB", "tcr") in got, "the beta-only record lost its TCR key"
    assert ("paired", "tcr") in got


def test_every_chain_metric_pair_spans_every_year(tmp_path, years):
    """R's `complete()` crosses the observed LEVELS; taking observed combinations drops pairs."""
    cum = panels.cumulative(_legacy(tmp_path), years)
    chains, metrics = cum["chains"].n_unique(), cum["metric"].n_unique()
    assert cum.height == chains * metrics * cum["year"].n_unique()


def test_totals_are_cumulative_and_never_decrease(tmp_path, years):
    cum = panels.cumulative(_legacy(tmp_path), years)
    for (_, _), g in cum.group_by(["chains", "metric"]):
        t = g.sort("year")["total"].to_list()
        assert t == sorted(t), "a cumulative count went down"


def test_by_year_draws_four_panels_with_the_callouts(tmp_path, years):
    cum = panels.cumulative(_legacy(tmp_path), years)
    ann = pl.DataFrame({"panel": ["tcr"], "year": [2021], "label": ["a mark"],
                        "hjust": [1], "vjust": [1]})
    fig = panels.by_year(cum, ann)
    assert len(fig.axes) == 4
    texts = [t.get_text() for ax in fig.axes for t in ax.texts]
    assert "a mark" in texts


def test_the_rcparams_are_the_manuscripts():
    # `2026-vdjdb-update` publishes with these; a dashboard panel and a paper panel should be the
    # same object. `fonttype 42` is what keeps PDF text editable rather than outlined.
    assert panels.RC["font.family"] == "Arial"
    assert panels.RC["pdf.fonttype"] == 42
    assert panels.RC["font.size"] == 7
