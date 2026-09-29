"""The profile has to add up. Nothing tested that, and it did not.

``motifs.tcrnet.background.*`` is timed inside ``motifs.tcrnet.enrichment``, so every row summed gave
a motif total of 268.7 s over a step that took 198 s and understated every ``share`` by 36 % - the
same ``share`` that ``tests/release/test_build_timings.py`` gates on, and the same network stages that
gate excludes from its denominator while the download also sat inside the enrichment row.
"""
from __future__ import annotations

import time

from vdjdb import timing


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
