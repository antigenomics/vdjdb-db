"""The build's own profile, gated on share of total rather than on seconds.

Marked ``release``: needs the timing report a real build writes (``VDJDB_TIMINGS``, default
``out/reports/build-timings.tsv``).

`ROADMAP_local.md` §49 profiled the assemble stage once, by hand, and found 156.19 s of 181.89 s in
`annotate.junction.add_junction_nt` - the measurement `antigenomics/vdjtools#181` rests on. Nothing in
the build recorded it, so the next stage to double would have been found the same way: by somebody
noticing the build felt slow.

**Why share and not seconds.** Seconds are a property of the host. A 4-vCPU runner is three to four
times slower than the laptop these numbers were first taken on, so a seconds bar is either useless or
fails on a busy runner - and a bar that fails on a correct build gets deleted, which is the lesson of
§52.2. Share is a ratio inside one run, so host speed cancels.

**What that cannot catch, stated rather than implied**: a *uniform* slowdown moves no share at all.
The absolute seconds are recorded in the artifact and printed into the step summary for exactly that
case, where a human comparing two runs is the instrument. What the gate catches is one stage blowing
up relative to the others, which is the failure that actually happened here (one call at 86 %).

The 10-point band is wide on purpose: `add_junction_nt` scales with cores while the polars stages do
not, so its share genuinely rises on a smaller host. 10 points absorbs that and still fails a stage
that goes from 1 % to 30 %.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.release

BASELINE = Path("rules/build_timings.tsv")


@pytest.fixture(scope="module")
def timings() -> pl.DataFrame:
    p = Path(os.environ.get("VDJDB_TIMINGS", "out/reports/build-timings.tsv"))
    if not p.exists():
        pytest.skip(f"no timing report at {p}; run `vdjdb build --out out/`")
    return pl.read_csv(p, separator="\t")


@pytest.fixture(scope="module")
def baseline() -> pl.DataFrame:
    if not BASELINE.exists():
        pytest.skip(f"no committed baseline at {BASELINE}")
    return pl.read_csv(BASELINE, separator="\t")


def test_every_stage_in_the_baseline_was_timed(timings, baseline) -> None:
    """A stage that stopped being timed is a stage nobody is watching, not a stage that got fast."""
    missing = set(baseline["stage"]) - set(timings["stage"])
    assert not missing, f"not timed by this build: {sorted(missing)}"


def test_no_stage_took_a_much_larger_share_than_its_baseline(timings, baseline) -> None:
    j = (baseline.join(timings.select("stage", pl.col("share").alias("now")), on="stage", how="left")
         .with_columns((pl.col("now") - pl.col("share")).alias("delta")))
    over = j.filter(pl.col("delta") > pl.col("tolerance"))
    assert over.height == 0, (
        "a stage grew by more than its tolerated share of the build; profile it, open an issue on the "
        f"repository whose code is slow, and re-record rules/build_timings.tsv:\n{over}")


def test_a_new_stage_is_recorded_before_it_can_dominate(timings, baseline) -> None:
    """An untimed stage taking a tenth of the build is the defect this file exists to prevent."""
    unknown = timings.filter(~pl.col("stage").is_in(baseline["stage"].to_list()))
    big = unknown.filter(pl.col("share") > 0.10)
    assert big.height == 0, (
        f"new stage(s) at over 10 % of the build and not in the baseline:\n{big}\n"
        "add them to rules/build_timings.tsv with their measured share")


def test_the_report_carries_the_input_size_and_the_core_count(timings) -> None:
    """A wall time without the input size and the core count is not a measurement (CLAUDE.md)."""
    assert (timings["rows"] > 0).all(), "the timing report must carry the record count"
    assert (timings["cores"] > 0).all(), "the timing report must carry the core count"
