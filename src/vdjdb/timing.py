"""Per-stage wall time, recorded by the build instead of profiled by hand afterwards.

The assemble stage's profile was measured once, on 2026-09-28, with an ad-hoc script: 156.19 s of
181.89 s in `annotate.junction.add_junction_nt`, 85.9 % of the whole (``ROADMAP_local.md`` §49). That
number is the reason `antigenomics/vdjtools#181` exists, and nothing in the build recorded it - so the
next time a stage doubles, someone has to notice the build feels slow and profile it again.

``CLAUDE.md`` asks for a profile as part of any bottleneck report, and asks for the wall time, the
share, the input size and the core count. This module makes all four fall out of an ordinary run.

**What is gated and what is only recorded, and why the difference is not laziness.** Absolute seconds
are a property of the host: a 4-vCPU runner is three to four times slower than the laptop these
numbers were first measured on, and gating them would either be so loose as to be useless or so tight
as to fail on a busy runner - which is how a timing bar gets deleted (``ROADMAP_local.md`` §52.2).
**Share of total** is a ratio measured inside one run, so host speed cancels, and it catches the thing
worth catching: one stage blowing up relative to the others. It cannot catch a uniform slowdown, and
that is stated rather than papered over.
"""
from __future__ import annotations

import time
from contextlib import contextmanager
from pathlib import Path

import polars as pl

#: ``stage -> seconds``, in call order. Module-level because the stages are spread over three
#: modules and threading a recorder through every signature would be a worse trade than one list.
#: Reset by :func:`reset`; a second build in the same process would otherwise append to the first.
_STAGES: list[tuple[str, float]] = []


def reset() -> None:
    _STAGES.clear()


@contextmanager
def stage(name: str):
    """Time one named stage. Nesting is allowed; the report says which are children."""
    start = time.perf_counter()
    try:
        yield
    finally:
        _STAGES.append((name, time.perf_counter() - start))


def frame(*, rows: int | None = None, cores: int | None = None) -> pl.DataFrame:
    """The timings as a table: one row per stage, with its share of the recorded total.

    ``share`` is of the summed leaf stages, not of the process, so it is comparable across hosts.
    ``rows`` and ``cores`` are carried because a wall time without the input size and the core count
    is not a measurement (``CLAUDE.md``).
    """
    import os

    total = sum(s for _, s in _STAGES) or 1.0
    return pl.DataFrame({
        "stage": [n for n, _ in _STAGES],
        "seconds": [round(s, 3) for _, s in _STAGES],
        "share": [round(s / total, 5) for _, s in _STAGES],
        "rows": [rows if rows is not None else -1] * len(_STAGES),
        "cores": [cores if cores is not None else (os.cpu_count() or -1)] * len(_STAGES),
    })


def write(path: Path, *, rows: int | None = None) -> pl.DataFrame:
    """Write the timing table beside the other reports. Returns it, for the caller to print."""
    df = frame(rows=rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path, separator="\t")
    return df


def report(df: pl.DataFrame) -> str:
    """The table as markdown, slowest first, for the CI step summary."""
    out = ["| stage | seconds | share | rows | cores |", "|---|--:|--:|--:|--:|"]
    for r in df.sort("seconds", descending=True).iter_rows(named=True):
        out.append(f"| `{r['stage']}` | {r['seconds']:.2f} | {r['share']:.1%} | "
                   f"{r['rows']:,} | {r['cores']} |")
    total = df["seconds"].sum()
    out.append(f"| **total** | **{total:.2f}** | | | |")
    return "\n".join(out)
