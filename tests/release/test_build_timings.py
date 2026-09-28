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

#: report -> (committed share baseline, peak-RSS budget in MiB). Both stages are covered, because
#: the pipeline's real memory peak is not in the one that was measured first.
REPORTS = {
    "build-timings.tsv": (Path("rules/build_timings.tsv"), 4096),
    "motif-timings.tsv": (Path("rules/motif_timings.tsv"), 10240),
}

#: Peak RSS budgets, MiB, per report. Measured 2026-09-28 on 16 cores:
#:
#: * **assemble 1,577 MiB**, flat across its stages because ``ru_maxrss`` is a high-water mark and the
#:   peak is set in the first one. Budget 4,096: 2.6x headroom.
#: * **motifs 6,898 MiB**, and this is the pipeline's real peak - 4.4x the assemble stage, reproducible
#:   at 6,871 and 6,898 over two runs, set by the two PWM-and-emit steps (1,005 -> 5,442 -> 6,898).
#:   Budget 10,240: 1.5x headroom, and 6 GB still free on a 16 GB runner.
#:
#: ⚠ **The first version of this file gated only the assemble report**, so it asserted 1,577 MiB
#: against a 16 GB runner while the pipeline actually peaked at 6,898 - 43 % of the runner rather than
#: the 10 % that figure implied (``ROADMAP_local.md`` section 58). Measuring one stage and calling it
#: the build is the mistake this parametrisation removes.
#:
#: Gated absolutely, unlike the seconds, because peak memory is a property of the data and the code
#: rather than of the host - the same build allocates the same way anywhere, and if anything a
#: *smaller* host allocates less, because polars chunks to fewer threads. So the laptop figure is the
#: conservative one to budget from.
#:
#: The claim these protect: the README put the pandas pipeline at 64 GB and the rewrite runs on a 16 GB
#: hosted runner. Until #610 **nothing checked it** - the `benchmark` mark promised "runtime and
#: peak-RSS budgets" and no test ever used it.


@pytest.fixture(params=sorted(REPORTS))
def report(request) -> tuple[str, pl.DataFrame, pl.DataFrame, int]:
    """``(name, measured, baseline, budget)`` for each timed stage of the pipeline."""
    name = request.param
    base_path, budget = REPORTS[name]
    got = Path(os.environ.get("VDJDB_REPORTS", "out/reports")) / name
    if not got.exists():
        pytest.skip(f"no {name}; run `vdjdb build` and `vdjdb motifs`")
    if not base_path.exists():
        pytest.skip(f"no committed baseline at {base_path}")
    return name, pl.read_csv(got, separator="\t"), pl.read_csv(base_path, separator="\t"), budget


def test_every_stage_in_the_baseline_was_timed(report) -> None:
    """A stage that stopped being timed is a stage nobody is watching, not a stage that got fast."""
    name, timings, baseline, _ = report
    missing = set(baseline["stage"]) - set(timings["stage"])
    assert not missing, f"{name}: not timed by this run: {sorted(missing)}"


def test_no_stage_took_a_much_larger_share_than_its_baseline(report) -> None:
    name, timings, baseline, _ = report
    j = (baseline.join(timings.select("stage", pl.col("share").alias("now")), on="stage", how="left")
         .with_columns((pl.col("now") - pl.col("share")).alias("delta")))
    over = j.filter(pl.col("delta") > pl.col("tolerance"))
    assert over.height == 0, (
        f"{name}: a stage grew by more than its tolerated share; profile it, open an issue on the "
        f"repository whose code is slow, and re-record its baseline:\n{over}")


def test_a_new_stage_is_recorded_before_it_can_dominate(report) -> None:
    """An untimed stage taking a tenth of the build is the defect this file exists to prevent."""
    name, timings, baseline, _ = report
    unknown = timings.filter(~pl.col("stage").is_in(baseline["stage"].to_list()))
    big = unknown.filter(pl.col("share") > 0.10)
    assert big.height == 0, (
        f"{name}: new stage(s) over 10 % of the run and not in the baseline:\n{big}\n"
        f"add them to {REPORTS[name][0]} with their measured share")


def test_the_build_stays_inside_its_memory_budget(report) -> None:
    """The claim is that this build runs on a 16 GB hosted runner. Now something checks it."""
    name, timings, _, budget = report
    peak = timings["peak_rss_mb"].max()
    assert peak > 0, f"{name} must carry a peak RSS; ru_maxrss is bytes on macOS and kilobytes on " \
                     f"Linux, and getting that wrong reads as a pass"
    assert peak <= budget, (
        f"{name}: peak RSS {peak:,.0f} MiB over the {budget:,} MiB budget. Profile what allocates, "
        f"and move the budget only with a measured reason.")


def test_the_report_carries_the_input_size_and_the_core_count(report) -> None:
    """A wall time without the input size and the core count is not a measurement (CLAUDE.md)."""
    _, timings, _, _ = report
    assert (timings["rows"] > 0).all(), "the timing report must carry the record count"
    assert (timings["cores"] > 0).all(), "the timing report must carry the core count"
