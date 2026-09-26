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


@pytest.mark.parametrize(("values", "expected"), [
    # `bw.nrd0(v)` in R 4.5.3, printed to 10 decimal places. Pinned against R's OUTPUT rather than
    # against the documented formula, because the two disagree: `bw.nrd` divides the IQR by 1.349
    # and `bw.nrd0` -- what `geom_density` actually defaults to -- divides by 1.34. Reading the
    # formula put the first implementation 0.67% out, which is close enough to look right.
    ([10, 11, 12, 12, 13, 13, 13, 14, 14, 15, 16, 18, 20], 1.2063415746),
    ([12] * 50 + [13] * 120 + [14] * 30, 0.0581931305),
])
def test_nrd0_matches_r(values, expected):
    assert panels.nrd0(values) == pytest.approx(expected, abs=1e-9)


def test_nrd0_falls_back_when_the_spread_is_zero():
    """R: `(lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)`. A constant vector must not give 0."""
    assert panels.nrd0([7.0] * 20) > 0
    assert panels.nrd0([0.0] * 20) > 0


def test_scores_draws_one_bar_group_per_class_and_chain(tmp_path):
    cohort = pl.DataFrame({
        "species": ["HomoSapiens"] * 6 + ["MusMusculus"],
        "mhc.class": ["MHCI"] * 3 + ["MHCII"] * 3 + ["MHCI"],
        "gene": ["TRA", "TRA", "TRB", "TRA", "TRB", "TRB", "TRA"],
        "vdjdb.score": ["0", "1", "0", "0", "2", "3", "0"],
    })
    fig = panels.scores(cohort)
    ax = fig.axes[0]
    assert [t.get_text() for t in ax.get_xticklabels()] == [
        "MHCI TRA", "MHCI TRB", "MHCII TRA", "MHCII TRB"]
    assert ax.get_yscale() == "log"


def test_v_hla_counts_source_rows_not_exploded_ones():
    """`length(unique(id))` in the R: the row id is assigned BEFORE the three explosions.

    One slim row can carry several MHC alleles and several V genes, so it lands in several cells --
    but it is one record in each, not one per combination. Counting after the explosion inflates
    every cell that a multi-allele row touches.
    """
    cohort = pl.DataFrame({
        "species": ["HomoSapiens"],
        "gene": ["TRB"],
        "mhc.class": ["MHCI"],
        "mhc.a": ["HLA-A*02:01,HLA-A*02:06"],   # two alleles that truncate to the SAME two-field
        "mhc.b": ["B2M"],
        "v.segm": ["TRBV7-9*01,TRBV7-9*03"],    # two alleles of the same V gene
    })
    cells = panels.v_hla(cohort, min_records=1)
    assert cells.height == 1, "the truncations collapse to one cell"
    assert cells["records"][0] == 1, "one source row must count once, not four times"
    assert cells["mhc"][0] == "A*02 / B2M"
    assert cells["v"][0] == "TRBV7-9"


def test_v_hla_drops_thin_alleles_and_genes():
    rows = {"species": [], "gene": [], "mhc.class": [], "mhc.a": [], "mhc.b": [], "v.segm": []}
    for n, (allele, v) in enumerate([("HLA-A*02:01", "TRBV1*01")] * 3
                                    + [("HLA-Z*99:01", "TRBV9*01")]):
        rows["species"].append("HomoSapiens")
        rows["gene"].append("TRB")
        rows["mhc.class"].append("MHCI")
        rows["mhc.a"].append(allele)
        rows["mhc.b"].append("B2M")
        rows["v.segm"].append(v)
        del n
    cells = panels.v_hla(pl.DataFrame(rows), min_records=3)
    assert cells["mhc"].to_list() == ["A*02 / B2M"]
    assert cells["records"][0] == 3


def test_epitope_length_draws_one_image_column_per_length():
    """Each bin is one image, not one artist per record -- the R stacks >100k bar segments."""
    cohort = pl.DataFrame({
        "species": ["HomoSapiens"] * 5,
        "mhc.class": ["MHCI"] * 3 + ["MHCII"] * 2,
        "cdr3": ["CASSA", "CASSBB", "CASSA", "CASSCCC", "CASSBB"],
        "antigen.epitope": ["GILGFVFTL", "GILGFVFTL", "NLVPMVATV",
                            "PKYVKQNTLKLAT", "PKYVKQNTLKLAT"],
    })
    fig = panels.epitope_length(cohort)
    assert len(fig.axes) == 2
    # one AxesImage per distinct epitope length present in that class
    assert len(fig.axes[0].images) == 1    # MHCI: only 9-mers
    assert len(fig.axes[1].images) == 1    # MHCII: only 13-mers
    assert [t.get_text() for t in fig.axes[0].get_xticklabels()] == ["9"]


def test_epitope_length_ignores_other_species():
    cohort = pl.DataFrame({
        "species": ["HomoSapiens", "MusMusculus"],
        "mhc.class": ["MHCI", "MHCI"],
        "cdr3": ["CASSA", "CASSB"],
        "antigen.epitope": ["GILGFVFTL", "SIINFEKL"],
    })
    fig = panels.epitope_length(cohort)
    assert [t.get_text() for t in fig.axes[0].get_xticklabels()] == ["9"]
