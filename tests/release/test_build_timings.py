"""Stage time budgets scale with the number of input rows on the same core count.

The reviewed reference profile records seconds and input rows. Assembly uses
records assembled from chunks; motifs use input chains. A stage's budget is
(reference seconds + tolerance * reference compute seconds) * current/reference rows * run scale.
The tolerance remains a fraction of reference compute time, and network fetches
contribute to neither the total nor the gated stages. This permits linear input
growth without editing the baseline. Multiply budgets by the median compute-stage
time-per-row ratio to adjust for common run speed. Uniform code and host slowdowns
cannot be distinguished; raw seconds remain recorded for review.
Core counts must match. Memory budgets remain absolute.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

pytestmark = pytest.mark.release

#: Fetch stages are excluded from both time budgets and the reference compute total.
NETWORK_STAGES = ("motifs.tcrnet.background.",)

#: Fixed four-core reference profiles from CI37919286987, before the gate migration.
#: Budgets retain the previous 0.1 allowances and absolute memory limits.
#: Refresh only for a reviewed algorithm or workload change, not ordinary row growth.
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
    """Compare stage budgets only on the recorded core count."""
    ran, recorded = int(timings["cores"][0]), int(baseline["cores"][0])
    if ran == recorded:
        return
    assert not os.environ.get("CI"), (
        f"{name}: the baseline was recorded on {recorded} cores and this CI run has {ran}. Time is "
        f"not comparable across host classes - use a reviewed matching reference profile.")
    pytest.skip(f"{name}: baseline recorded on {recorded} cores, this run has {ran}; the timing "
                f"comparison needs one host class. The memory budget still applies.")


def _compute_only(df: pl.DataFrame) -> pl.DataFrame:
    """Compute-only shares identify previously unrecorded dominant stages."""
    d = df.filter(~pl.col("stage").str.starts_with(NETWORK_STAGES[0]))
    for prefix in NETWORK_STAGES[1:]:
        d = d.filter(~pl.col("stage").str.starts_with(prefix))
    return d.with_columns(pl.col("share") / pl.col("share").sum())


def test_no_stage_exceeds_its_row_scaled_time_budget(report) -> None:
    from vdjdb.timing import row_scaled_limits

    name, timings, baseline, _ = report
    _same_host(name, timings, baseline)
    limits = row_scaled_limits(timings, baseline, exclude_prefixes=NETWORK_STAGES)
    over = limits.filter(pl.col("now_seconds") > pl.col("budget_seconds"))
    assert over.height == 0, (
        f"{name}: stage time exceeded its input-row-scaled budget. Profile it and "
        f"open an issue in the repository owning the slow code:\n{over}")


def test_a_new_stage_is_recorded_before_it_can_dominate(report) -> None:
    """An untimed stage taking a tenth of the build is the defect this file exists to prevent.

    Compute-only for the same reason the share comparison is: `background.HomoSapiens.TRB` is 0.209
    of the raw total, so a *new* background - one more species or chain in `chunks/` - would clear
    0.10 on download time alone and red the build for a fetch nobody wrote.
    """
    name, timings, baseline, _ = report
    timings = _compute_only(timings)
    unknown = timings.filter(~pl.col("stage").is_in(baseline["stage"].to_list()))
    big = unknown.filter(pl.col("share") > 0.10)
    assert big.height == 0, (
        f"{name}: new stage(s) over 10 % of the run and not in the baseline:\n{big}\n"
        f"add them to {REPORTS[name][0]} with their measured seconds and input rows")


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


def test_assessment_process_tree_stays_inside_its_memory_budget(report) -> None:
    name, timings, _, _ = report
    if name != "build-timings.tsv":
        return
    assessment = timings.filter(pl.col("stage") == "build_epitope_assessment")
    assert assessment.height == 1
    assert "peak_tree_rss_mb" in assessment.columns
    peak = assessment["peak_tree_rss_mb"].cast(pl.Float64).item()
    # Cold four-process assessment measured 4,931.4 MiB across parent and workers.
    # Keep the existing 4,096 MiB parent gate and separately budget the process tree.
    assert 0 < peak <= 8192, f"assessment process-tree peak {peak:,.1f} MiB exceeds 8,192 MiB"
