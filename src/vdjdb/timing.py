"""Exclusive per-stage wall time and process memory measured during each build.

Reports retain seconds, shares, input rows and core count. Fixed reference
profiles define time budgets scaled by input rows on the same core count.
Network fetches are excluded from compute budgets; memory limits are absolute.
"""
from __future__ import annotations

import os
import resource
import sys
import threading
import time
from contextlib import contextmanager, suppress
from pathlib import Path

import polars as pl
import psutil


def row_scaled_limits(timings: pl.DataFrame, baseline: pl.DataFrame, *,
                      exclude_prefixes: tuple[str, ...] = ()) -> pl.DataFrame:
    """Compare stage seconds with budgets scaled by the recorded input row count.

    A tolerance is a fraction of the reference's total compute time, preserving
    the former percentage-point allowance. Fetch time contributes to neither
    the reference total nor the gated stages. Row counts describe records for
    assembly and chains for motifs; they must agree within each profile. The
    median compute-stage time-per-row ratio adjusts for common run speed.
    Uniform code and host slowdowns cannot be distinguished by this profile;
    the raw seconds remain available for that review.
    """
    for name, frame in (("timings", timings), ("baseline", baseline)):
        required = {"stage", "seconds", "rows", "cores"}
        if name == "baseline":
            required.add("tolerance")
        if not required <= set(frame.columns):
            raise ValueError(f"{name}: missing timing budget columns")
        if (frame.is_empty() or frame.select(sorted(required)).null_count().sum_horizontal().item()
                or frame["stage"].n_unique() != frame.height
                or frame["rows"].n_unique() != 1 or frame["cores"].n_unique() != 1
                or not (frame["rows"] > 0).all() or not (frame["cores"] > 0).all()
                or not frame["seconds"].is_finite().all()
                or not (frame["seconds"] >= 0).all()):
            raise ValueError(f"{name}: invalid timing profile")
    if timings["cores"][0] != baseline["cores"][0]:
        raise ValueError("timing budgets require the same core count")
    if (not baseline["tolerance"].is_finite().all()
            or not (baseline["tolerance"] >= 0).all()):
        raise ValueError("baseline: invalid timing tolerance")
    reference, current = baseline, timings
    for prefix in exclude_prefixes:
        reference = reference.filter(~pl.col("stage").str.starts_with(prefix))
        current = current.filter(~pl.col("stage").str.starts_with(prefix))
    total = reference["seconds"].sum()
    if not total or total <= 0:
        raise ValueError("baseline: no positive compute time")
    if not set(reference["stage"]) <= set(current["stage"]):
        raise ValueError("timings: a reference stage was not measured")
    scale = timings["rows"][0] / baseline["rows"][0]
    compared = reference.select("stage", "seconds", "tolerance").join(
        current.select("stage", pl.col("seconds").alias("now_seconds")),
        on="stage", how="left")
    speed = (compared.filter(pl.col("seconds") > 0)
             .select((pl.col("now_seconds") / (pl.col("seconds") * scale)).median())
             .item())
    if not speed or speed <= 0:
        raise ValueError("timings: no positive common run speed")
    return (compared
            .with_columns(
                (pl.col("seconds") * scale).alias("expected_seconds"),
                ((pl.col("seconds") + pl.col("tolerance") * total) * scale * speed)
                .alias("budget_seconds"),
                pl.lit(scale).alias("row_scale"),
                pl.lit(speed).alias("run_scale"))
            .sort("stage"))


def peak_rss_mb() -> float:
    """The process high-water mark in MiB.

    ``ru_maxrss`` is **kilobytes on Linux and bytes on macOS** - the one portability trap here, and
    getting it wrong would report the CI runner using a thousandth of its real memory, which is the
    direction that reads as a pass.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw / (1024 * 1024) if sys.platform == "darwin" else raw / 1024


#: One entry per stage, in the order the stages were *entered*: ``[name, parent index, seconds spent
#: inside it, seconds spent inside its timed children, peak RSS in MiB on exit]``. Module-level
#: because the stages are spread over three modules and threading a recorder through every signature
#: would be a worse trade than one list.
#: Reset by :func:`reset`; a second build in the same process would otherwise append to the first.
_STAGES: list[list] = []

#: Indices into :data:`_STAGES` of the stages currently open, outermost first. This is what makes a
#: nested stage attributable: without it a child's seconds are counted twice, once in its own row and
#: once inside its parent's.
_OPEN: list[int] = []


def reset() -> None:
    _STAGES.clear()
    _OPEN.clear()


@contextmanager
def stage(name: str, *, process_tree: bool = False):
    """Time one named stage, exclusive of any stage timed inside it.

    **A nested stage used to be counted twice.** ``motifs.tcrnet.background.*`` runs inside
    ``motifs.tcrnet.enrichment``, so summing every row gave a motif total of 268.7 s over a step that
    took 198 s, and every ``share`` in the report was understated by 36 % - including the shares
    ``tests/release/test_build_timings.py`` gates on, and including the network shares that gate
    deliberately excludes from its denominator, which it could not actually exclude while the
    download was also inside the enrichment row. Each stage now reports the time spent in it and
    *not* in a timed child, so the rows sum to the wall clock exactly once and ``parent`` says where
    a child belongs.
    """
    me = len(_STAGES)
    _STAGES.append([name, _OPEN[-1] if _OPEN else -1, 0.0, 0.0, 0.0, None])
    _OPEN.append(me)
    stop = threading.Event()
    peak = [0]
    errors = []

    def sample():
        try:
            parent = psutil.Process()
            processes = [parent, *parent.children(recursive=True)]
            rss = 0
            for process in processes:
                # A child can exit between enumeration and reading its RSS.
                with suppress(psutil.NoSuchProcess):
                    rss += process.memory_info().rss
            peak[0] = max(peak[0], rss)
        except Exception as error:
            errors.append(error)
            stop.set()

    def monitor():
        while not stop.wait(0.05):
            sample()

    sampler = None
    if process_tree:
        sample()
        sampler = threading.Thread(target=monitor, name="vdjdb-rss", daemon=True)
        sampler.start()
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        if sampler is not None:
            stop.set()
            sampler.join()
            sample()
            _STAGES[me][5] = peak[0] / (1024 * 1024)
        _OPEN.pop()
        # ``ru_maxrss`` is a high-water mark, so the value after a stage is "peak so far" rather than
        # that stage's own footprint. That is the useful reading: it says which stage raised the peak.
        _STAGES[me][2], _STAGES[me][4] = elapsed, peak_rss_mb()
        if _STAGES[me][1] >= 0:
            _STAGES[_STAGES[me][1]][3] += elapsed
        if errors:
            raise RuntimeError("process-tree RSS measurement failed") from errors[0]


def cores_available() -> int:
    """CPUs this process may actually use.

    ``os.cpu_count()`` reports the machine, not the allocation, so under a SLURM ``-c 4`` cpuset or a
    container quota it says 40 where four are usable - and the share baseline is keyed on this number
    (``tests/release/test_build_timings.py``), so an overstated count makes the gate refuse to compare
    on a host that is in fact the same class as the runner. ``sched_getaffinity`` is Linux-only;
    macOS has no cpuset to misreport.
    """
    getaffinity = getattr(os, "sched_getaffinity", None)
    return len(getaffinity(0)) if getaffinity else (os.cpu_count() or -1)


def _self_seconds() -> list[float]:
    """Per stage, the seconds spent in it and not in a stage timed inside it."""
    return [max(e[2] - e[3], 0.0) for e in _STAGES]


def frame(*, rows: int | None = None, cores: int | None = None) -> pl.DataFrame:
    """The timings as a table: one row per stage, with its share of the recorded total.

    ``seconds`` is exclusive of timed children and ``parent`` names the stage a row sits inside, so
    the column sums to the wall clock and ``share`` sums to 1. A stage's inclusive duration is its
    own seconds plus its children's, which is why the parent is recorded rather than the total.

    ``rows`` and ``cores`` are carried because a wall time without the input size and the core count
    is not a measurement (``CLAUDE.md``).
    """
    self_s = _self_seconds()
    total = sum(self_s) or 1.0
    return pl.DataFrame({
        "stage": [e[0] for e in _STAGES],
        "parent": [_STAGES[e[1]][0] if e[1] >= 0 else "" for e in _STAGES],
        "seconds": [round(s, 3) for s in self_s],
        "share": [round(s / total, 5) for s in self_s],
        "peak_rss_mb": [round(e[4], 1) for e in _STAGES],
        "peak_tree_rss_mb": [f"{e[5]:.1f}" if e[5] is not None else "" for e in _STAGES],
        "rows": [rows if rows is not None else -1] * len(_STAGES),
        "cores": [cores if cores is not None else cores_available()] * len(_STAGES),
    })


def write(path: Path, *, rows: int | None = None) -> pl.DataFrame:
    """Write the timing table beside the other reports. Returns it, for the caller to print."""
    df = frame(rows=rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path, separator="\t")
    return df


def report(df: pl.DataFrame) -> str:
    """The table as markdown, slowest first, for the CI step summary.

    ``inside`` is the parent stage, empty for a top-level one. Seconds are exclusive of children, so
    the total is the wall clock and not the 36 % overstatement summing every row used to give.
    """
    out = ["| stage | inside | seconds | share | parent peak RSS MiB | "
           "sampled tree peak RSS MiB | rows | cores |",
           "|---|---|--:|--:|--:|--:|--:|--:|"]
    for r in df.sort("seconds", descending=True).iter_rows(named=True):
        parent = r.get("parent") or ""
        tree = r.get("peak_tree_rss_mb")
        tree_text = f"{float(tree):,.0f}" if tree else "-"
        out.append(f"| `{r['stage']}` | {f'`{parent}`' if parent else ''} | {r['seconds']:.2f} | "
                   f"{r['share']:.1%} | {r['peak_rss_mb']:,.0f} | {tree_text} | "
                   f"{r['rows']:,} | {r['cores']} |")
    trees = [float(v) for v in df['peak_tree_rss_mb'] if v]
    tree_total = f"{max(trees):,.0f}" if trees else "-"
    out.append(f"| **total** | | **{df['seconds'].sum():.2f}** | | "
               f"**{df['peak_rss_mb'].max():,.0f}** | **{tree_total}** | | |")
    return "\n".join(out)
