"""The profile has to add up. Nothing tested that, and it did not.

``motifs.tcrnet.background.*`` is timed inside ``motifs.tcrnet.enrichment``, so every row summed gave
a motif total of 268.7 s over a step that took 198 s and understated every ``share`` by 36 % - the
same ``share`` that ``tests/release/test_build_timings.py`` gates on, and the same network stages that
gate excludes from its denominator while the download also sat inside the enrichment row.
"""
from __future__ import annotations

import subprocess
import sys
import threading
import time

import polars as pl
import psutil
import pytest

from vdjdb import timing


def _profile(rows=100, *, seconds=(80.0, 20.0, 10.0), network=5.0, cores=4):
    return pl.DataFrame({
        "stage": ["score", "emit", "sort", "fetch.background"],
        "seconds": [*seconds, network], "rows": [rows] * 4,
        "cores": [cores] * 4, "tolerance": [0.1] * 4,
    })


def test_timing_budgets_scale_with_input_growth_and_shrinkage():
    for rows, seconds in ((200, (160.0, 40.0, 20.0)), (50, (40.0, 10.0, 5.0))):
        limits = timing.row_scaled_limits(_profile(rows, seconds=seconds), _profile(),
                                          exclude_prefixes=("fetch.",))
        assert limits["budget_seconds"].to_list() == [31 * rows / 100, 91 * rows / 100, 21 * rows / 100]
        assert not limits.filter(pl.col("now_seconds") > pl.col("budget_seconds")).height


def test_row_scaling_catches_individual_stage_regressions():
    limits = timing.row_scaled_limits(_profile(200, seconds=(190.0, 40.0, 20.0)),
                                      _profile(), exclude_prefixes=("fetch.",))
    assert limits.filter(pl.col("now_seconds") > pl.col("budget_seconds"))["stage"].to_list() == ["score"]


def test_common_speed_variation_is_recorded_without_moving_the_reference():
    limits = timing.row_scaled_limits(_profile(200, seconds=(192.0, 48.0, 24.0)),
                                      _profile(), exclude_prefixes=("fetch.",))
    assert limits["row_scale"].to_list() == [2.0] * 3
    assert limits["run_scale"].to_list() == pytest.approx([1.2] * 3)
    assert not limits.filter(pl.col("now_seconds") > pl.col("budget_seconds")).height


def test_fetch_seconds_do_not_change_compute_budgets():
    limits = timing.row_scaled_limits(_profile(network=50000.0), _profile(network=70000.0),
                                      exclude_prefixes=("fetch.",))
    assert limits["stage"].to_list() == ["emit", "score", "sort"]
    assert limits["budget_seconds"].to_list() == [31.0, 91.0, 21.0]


@pytest.mark.parametrize("bad", [
    _profile(rows=0), _profile(cores=2), _profile(seconds=(float("nan"), 20.0, 10.0)),
    _profile().with_columns(pl.Series("rows", [100, 200, 100, 100])),
    _profile().with_columns(pl.lit(None, dtype=pl.Int64).alias("rows")),
    _profile().with_columns(pl.lit(-0.1).alias("tolerance")),
])
def test_invalid_timing_reference_is_rejected(bad):
    with pytest.raises(ValueError):
        timing.row_scaled_limits(_profile(), bad, exclude_prefixes=("fetch.",))


def _sleepy(seconds: float) -> None:
    end = time.perf_counter() + seconds
    while time.perf_counter() < end:
        pass


def test_a_nested_stage_is_not_counted_in_its_parent() -> None:
    timing.reset()
    with timing.stage("parent"):
        _sleepy(0.02)
        with timing.stage("child"):
            _sleepy(0.05)
    df = timing.frame(rows=1).sort("stage")
    assert df["stage"].to_list() == ["child", "parent"]
    assert df["parent"].to_list() == ["parent", ""]
    # The child is the larger of the two, which is the whole reading the old report inverted.
    assert df.filter(stage="child")["seconds"][0] > df.filter(stage="parent")["seconds"][0]
    # `share` is rounded to five places, so the sum is exact only to that.
    assert abs(df["share"].sum() - 1.0) < 1e-4


def test_the_total_is_the_wall_clock() -> None:
    timing.reset()
    start = time.perf_counter()
    with timing.stage("outer"):
        with timing.stage("inner"):
            _sleepy(0.03)
        _sleepy(0.01)
    wall = time.perf_counter() - start
    total = timing.frame(rows=1)["seconds"].sum()
    assert abs(total - wall) < 0.02, f"{total:.3f} s recorded against {wall:.3f} s of wall clock"


def test_sibling_stages_still_sum_to_one() -> None:
    timing.reset()
    for name in ("a", "b", "c"):
        with timing.stage(name):
            _sleepy(0.01)
    df = timing.frame(rows=1)
    assert df["parent"].to_list() == ["", "", ""]
    # `share` is rounded to five places, so the sum is exact only to that.
    assert abs(df["share"].sum() - 1.0) < 1e-4


def test_reset_clears_an_interrupted_stack() -> None:
    """A raising stage leaves the stack unwound; a second build must not inherit a parent."""
    timing.reset()
    try:
        with timing.stage("boom"):
            raise RuntimeError
    except RuntimeError:
        pass
    timing.reset()
    with timing.stage("fresh"):
        _sleepy(0.001)
    assert timing.frame(rows=1)["parent"].to_list() == [""]


def test_process_tree_sampling_includes_a_live_worker():
    timing.reset()
    before = psutil.Process().memory_info().rss / (1024 * 1024)
    with timing.stage('workers', process_tree=True), subprocess.Popen([sys.executable, '-c',
                           'import time; data=bytearray(128*1024*1024); '
                           'print("ready",flush=True); time.sleep(.3)'],
                          stdout=subprocess.PIPE, text=True) as child:
        assert child.stdout.readline().strip() == 'ready'
        assert child.wait(timeout=5) == 0
    peak = float(timing.frame()['peak_tree_rss_mb'].item())
    assert peak > before + 96
    assert not any(t.name == 'vdjdb-rss' for t in threading.enumerate())


def test_process_tree_measurement_failure_is_explicit_and_cleans_up(monkeypatch):
    timing.reset()

    def denied():
        raise psutil.AccessDenied()

    monkeypatch.setattr(psutil, 'Process', denied)
    with pytest.raises(RuntimeError, match='process-tree RSS measurement failed'), \
            timing.stage('denied', process_tree=True):
        pass
    assert not timing._OPEN
    assert not any(t.name == 'vdjdb-rss' for t in threading.enumerate())
