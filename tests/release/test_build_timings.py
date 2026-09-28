"""The build's own profile, gated on share of total rather than on seconds.

Marked ``release``: needs the timing report a real build writes (``VDJDB_TIMINGS``, default
``out/reports/build-timings.tsv``).

`ROADMAP_local.md` §49 profiled the assemble stage once, by hand, and found 156.19 s of 181.89 s in
`annotate.junction.add_junction_nt` - the measurement `antigenomics/vdjtools#181` rests on. Nothing in
the build recorded it, so the next stage to double would have been found the same way: by somebody
noticing the build felt slow.

**Why share and not seconds.** Seconds are a property of the host. A 4-vCPU runner is two to five
times slower than the laptop these numbers were first taken on, so a seconds bar is either useless or
fails on a busy runner - and a bar that fails on a correct build gets deleted, which is the lesson of
§52.2. Share is a ratio inside one run, so a *uniform* host slowdown cancels.

**Share cancels host speed only when the stages scale alike, and that was measured on both reports
rather than assumed.** Same corpus, 16-core laptop against the 4-vCPU runner:

* **assemble**: 187.3 s against 443.2 s, 2.37x, and the largest share difference over nine stages is
  **1.3 points**. Every stage there is slower by about the same factor, so the ratio does cancel.
* **motifs**: 40.1 s against 211.3 s, 5.27x, and the largest share difference over thirteen stages is
  **33.8 points** - because the per-stage slowdown ranges from 1.46x (`tcrnet.pwm_and_emit`, polars)
  to 24.0x (`tcrnet.background`, which streams four backgrounds from HuggingFace on the runner and
  reads them out of the local cache on a laptop).

So the share comparison runs **only when the run's core count matches the count the baseline was
recorded on**, which the baseline now carries. Under CI that is an assertion rather than a skip: the
runner is fixed at 4 vCPU, so a mismatch means the baseline was recorded somewhere else and has to be
re-recorded, and a gate that skips itself is the failure mode that reads as a pass (§59.2).

**What share cannot catch, stated rather than implied**: a uniform slowdown moves no share at all. The
absolute seconds are recorded in the artifact and printed into the step summary for exactly that case,
where a human comparing two runs is the instrument. What the gate catches is one stage blowing up
relative to the others, which is the failure that actually happened here (one call at 86 %).

**Peak RSS needs none of this**, so it is gated absolutely and on every host. Five measurements of
the motif stage: 6,871 and 6,898 MiB on the laptop, 6,786, 6,802 and 7,531 MiB on the runner. The
runner spread is 11 %, so it is not the invariant the first pair suggested, and the budget is set
against the largest observation rather than the flattering one.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.release

#: Stages that are a network fetch rather than a computation. They are timed and recorded - the
#: runner spends 55.6 s and 73.0 s of the motif stage acquiring four backgrounds over two measured
#: runs - but they are **excluded from the share comparison, denominator included**, because their
#: duration is GitHub's network and not this repository's code.
#:
#: Without that, the gate is one slow download from red. Measured against the committed baseline,
#: where the fetch is 0.26803 of the recorded total: if the fetch were free, `tcrnet.enrichment`
#: renormalises 0.29646 -> 0.40502, delta +0.10856; if it took twice as long,
#: `background.HomoSapiens.TRB` goes 0.20905 -> 0.32972, delta +0.12067. Both exceed the 0.1 band
#: while nothing about the code changed, and the fetch already varied 55.6 s -> 73.0 s (1.31x)
#: between two runs of one commit.
#:
#: Compared compute-only, the same two runs differ by at most **0.05290** and the fetch moves nothing.
NETWORK_STAGES = ("motifs.tcrnet.background.",)

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
#: * **motifs 6,786 to 7,531 MiB**, and this is the pipeline's peak - 4.4x the assemble stage, set by
#:   the two PWM-and-emit steps (1,005 -> 5,442 -> 6,898 on the laptop). Five measurements: 6,871 and
#:   6,898 on a 16-core laptop, 6,786, 6,802 and 7,531 on the runner. Budget 10,240: **1.36x headroom
#:   against the largest**, and 2.6 GB still free on a 16 GB runner.
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
        # Under CI a missing report is the defect, not a reason to stand down: the build ran, so the
        # report exists or the stage stopped writing it. `vdjdb motifs` wrote this one to
        # `<--out>/reports/` for one commit, which put it at `out/motifs/reports/` on the runner and
        # `out/reports/` on a laptop, and a skip would have made the memory gate silently absent from
        # every CI run while passing locally.
        assert not os.environ.get("CI"), (
            f"no {name} at {got} after a CI build; the stage that writes it either did not run or "
            f"writes somewhere else now")
        pytest.skip(f"no {name}; run `vdjdb build` and `vdjdb motifs`")
    if not base_path.exists():
        pytest.skip(f"no committed baseline at {base_path}")
    return name, pl.read_csv(got, separator="\t"), pl.read_csv(base_path, separator="\t"), budget


def test_every_stage_in_the_baseline_was_timed(report) -> None:
    """A stage that stopped being timed is a stage nobody is watching, not a stage that got fast."""
    name, timings, baseline, _ = report
    missing = set(baseline["stage"]) - set(timings["stage"])
    assert not missing, f"{name}: not timed by this run: {sorted(missing)}"


def _same_host(name: str, timings: pl.DataFrame, baseline: pl.DataFrame) -> None:
    """Share is comparable within a host class, not across them. See the module docstring."""
    ran, recorded = int(timings["cores"][0]), int(baseline["cores"][0])
    if ran == recorded:
        return
    assert not os.environ.get("CI"), (
        f"{name}: the baseline was recorded on {recorded} cores and this CI run has {ran}. Share is "
        f"not comparable across host classes - re-record the baseline from a CI run.")
    pytest.skip(f"{name}: baseline recorded on {recorded} cores, this run has {ran}; the share "
                f"comparison needs one host class. The memory budget still applies.")


def _compute_only(df: pl.DataFrame) -> pl.DataFrame:
    """Drop the network stages and renormalise ``share`` over what is left. See NETWORK_STAGES.

    Both sides of the comparison go through this, so the committed baseline keeps the shares as they
    were measured and the fetch is excluded from the denominator on the run as well.
    """
    d = df.filter(~pl.col("stage").str.starts_with(NETWORK_STAGES[0]))
    for prefix in NETWORK_STAGES[1:]:
        d = d.filter(~pl.col("stage").str.starts_with(prefix))
    return d.with_columns(pl.col("share") / pl.col("share").sum())


def test_no_stage_took_a_much_larger_share_than_its_baseline(report) -> None:
    name, timings, baseline, _ = report
    _same_host(name, timings, baseline)
    timings, baseline = _compute_only(timings), _compute_only(baseline)
    j = (baseline.join(timings.select("stage", pl.col("share").alias("now")), on="stage", how="left")
         .with_columns((pl.col("now") - pl.col("share")).alias("delta")))
    over = j.filter(pl.col("delta") > pl.col("tolerance"))
    assert over.height == 0, (
        f"{name}: a stage grew by more than its tolerated share; profile it, open an issue on the "
        f"repository whose code is slow, and re-record its baseline:\n{over}")


def test_a_new_stage_is_recorded_before_it_can_dominate(report) -> None:
    """An untimed stage taking a tenth of the build is the defect this file exists to prevent.

    Compute-only for the same reason the share comparison is: `background.HomoSapiens.TRB` is 0.209
    of the raw total, so a *new* background - one more species or chain in `chunks/` - would clear
    0.10 on download time alone and red the build for a fetch nobody wrote.
    """
    name, timings, baseline, _ = report
    timings, baseline = _compute_only(timings), _compute_only(baseline)
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
