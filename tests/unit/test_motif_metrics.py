"""The motif-metric gate: does it actually catch a regression, in both directions?

A gate that cannot fail is worse than no gate, because it reads as a passing check. Every test here
builds a metric frame by hand and asserts the gate's verdict on it, so none of them needs a build.

The committed baseline itself is checked too: it has to cover the three sources the report promises
(`docs/clustering.md`), or a source could stop being produced and the missing-row rule would be the
only thing left to notice.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.validate import motif_metrics as mm

#: The frame shape every gate joins on. An empty frame built without it carries no columns at all,
#: which fails the join rather than reporting the missing source.
SCHEMA = mm.COLUMNS


def _rows(**by_axis: float) -> list[dict]:
    return [{"species": "HomoSapiens", "gene": "TRB", "source": src, "axis": axis,
             "direction": mm.AXES[axis][0], "tolerance": mm.AXES[axis][1], "value": v,
             "reason": None}
            for src, vals in by_axis.items() for axis, v in vals.items()]  # type: ignore[union-attr]


def frame(**by_source: dict) -> pl.DataFrame:
    return pl.DataFrame(_rows(**by_source), schema=SCHEMA)


def test_a_fall_in_an_up_axis_beyond_tolerance_is_a_regression() -> None:
    m = frame(latest={"purity": 0.9790}, **{"current-tcrnet": {"purity": 0.9600}})
    bad = mm.regressions_against_latest(m)
    assert bad.height == 1
    assert bad["axis"][0] == "purity"
    assert bad["delta"][0] == pytest.approx(-0.019, abs=1e-9)


def test_a_fall_inside_tolerance_is_not_a_regression() -> None:
    """0.005 is room for the corpus to grow between a recorded value and a run, not for noise."""
    m = frame(latest={"purity": 0.9790}, **{"current-tcrnet": {"purity": 0.9750}})
    assert mm.regressions_against_latest(m).height == 0


def test_percolation_is_recorded_and_never_gated() -> None:
    """Percolation is an epitope's largest-cluster share, and what is right depends on the epitope:
    one prominent motif (A*02 GIL) should percolate, a featureless one (A*02 NLV) should not. So
    neither a rise nor a fall is a regression; fusion of different motifs shows in purity."""
    assert mm.AXES["percolation_median"][0] == "record"
    assert mm.AXES["percolation_excess"][0] == "record"
    higher = frame(latest={"percolation_median": 0.6667, "percolation_excess": 0.0},
                   **{"current-tcrnet": {"percolation_median": 0.9000, "percolation_excess": 0.2333}})
    lower = frame(latest={"percolation_median": 0.6667},
                  **{"current-tcrnet": {"percolation_median": 0.1000}})
    assert mm.regressions_against_latest(higher).height == 0
    assert mm.regressions_against_latest(lower).height == 0


def test_losing_one_covered_epitope_is_a_regression() -> None:
    """`epitopes` has tolerance 0: it counts epitopes that get any denoising at all."""
    m = frame(latest={"epitopes": 109.0}, **{"current-tcrnet": {"epitopes": 108.0}})
    assert mm.regressions_against_latest(m).height == 1


def test_a_recorded_axis_never_gates() -> None:
    """More clusters is better for coverage and worse for parsimony, so neither direction is a
    regression on its own."""
    m = frame(latest={"clusters": 1051.0}, **{"current-tcrnet": {"clusters": 10.0}})
    assert mm.regressions_against_latest(m).height == 0


def test_baseline_drift_is_reported_separately_from_the_build_gate() -> None:
    """The legacy release's own numbers move when the cohort changes, and that is the case the
    2026-09-28 chunk merges hit. It is a drift against the baseline, not a build regression."""
    measured = frame(legacy={"retention": 0.3159}, latest={"retention": 0.3159})
    base = frame(legacy={"retention": 0.3218}, latest={"retention": 0.3218})
    assert mm.regressions_against_latest(measured).height == 0
    drift = mm.regressions_against_baseline(measured, base)
    # Both released rows move together, because they are the same fixed file scored on a cohort that
    # changed. Reporting only one of them would be the misleading half of the story.
    assert sorted(drift["source"]) == ["latest", "legacy"]
    assert drift["delta"].to_list() == [pytest.approx(-0.0059, abs=1e-9)] * 2


def test_a_source_that_stopped_being_produced_fails_rather_than_passing() -> None:
    """A left join would leave a vanished clustering scoring nothing, which reads as no regression."""
    base = frame(**{"current-tcremp": {"purity": 0.9829}})
    assert mm.regressions_against_baseline(frame(), base).height == 1


def test_the_report_names_every_axis_and_both_comparisons() -> None:
    m = frame(latest={"purity": 0.9790}, **{"current-tcrnet": {"purity": 0.9800}})
    text = mm.report(m, frame(latest={"purity": 0.9790}))
    assert "`purity`" in text and "current-tcrnet" in text
    assert "+0.0010" in text          # vs latest
    assert text.count("\n") == 3      # header, rule, two rows


@pytest.mark.skipif(not mm.BASELINE.exists(), reason="baseline not recorded yet")
def test_the_committed_baseline_covers_every_source_and_gated_axis() -> None:
    b = mm.load_baseline()
    assert set(b["source"].unique()) == set(mm.SOURCES), "a source is missing from the baseline"
    for (gene, source), grp in b.group_by(["gene", "source"]):
        missing = set(mm.AXES) - set(grp["axis"])
        assert not missing, f"{gene} {source} has no baseline for {sorted(missing)}"
    assert b.filter(pl.col("direction") == "up").height, "nothing is gated upward"


def test_a_declared_trade_is_accepted_and_a_worse_one_is_not() -> None:
    """An axis that is worse on purpose is declared with its reason and capped at the declared value.
    The same axis drifting further is not accepted, and must still fail.
    """
    declared = frame(**{"current-tcremp": {"retention": 0.2000}}).with_columns(
        pl.lit("buys a wider ball and two epitopes").alias("reason"))
    at = frame(latest={"retention": 0.2500}, **{"current-tcremp": {"retention": 0.2000}})
    worse = frame(latest={"retention": 0.2500}, **{"current-tcremp": {"retention": 0.1500}})
    assert mm.regressions_against_latest(at, declared).height == 0
    assert mm.regressions_against_latest(at).height == 1, "undeclared, it must be reported"
    assert mm.regressions_against_latest(worse, declared).height == 1
