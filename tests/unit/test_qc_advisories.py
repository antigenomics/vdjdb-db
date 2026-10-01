"""The QC advisory counts, recorded and gated. `chunks/` is the data, so these are its quality signal.

`vdjdb qc` prints its per-rule counts and, until now, wrote them nowhere: an advisory that jumped from
209 rows to 2,000 on a chunk merge passed every gate and appeared only in a log line nobody diffs.
`--strict` fails on a *fatal* finding, which is the right behaviour and a different question - an
advisory is advisory because only a curator can resolve it, not because its size does not matter.

The repo's own rule already asks a chunk commit to state "what the build shows". This makes that
mechanical: the counts are an artifact, and a rise beyond its declared tolerance fails until somebody
moves `rules/qc_advisories.tsv` in a commit that says why.

**Tolerances differ by rule on purpose.** `duplicate` counts within-chunk duplicate rows and grows with
the corpus, so its band is wide; a malformed column name is a one-line repair in the submission and its
band is zero. A single number for both would either wave through a real jump or fail every curation PR.

Unit-marked, not release-marked: it needs `chunks/` and the QC report, which every chunk pull request
already produces, and this is the tier that runs there.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

BASELINE = Path("rules/qc_advisories.tsv")


@pytest.fixture(scope="module")
def summary() -> pl.DataFrame:
    p = Path(os.environ.get("VDJDB_QC_SUMMARY", "out/reports/qc-summary.tsv"))
    if not p.exists():
        pytest.skip(f"no QC summary at {p}; run `vdjdb qc --report out/reports/qc.tsv`")
    return pl.read_csv(p, separator="\t")


@pytest.fixture(scope="module")
def baseline() -> pl.DataFrame:
    if not BASELINE.exists():
        pytest.skip(f"no committed baseline at {BASELINE}")
    return pl.read_csv(BASELINE, separator="\t")


def test_no_advisory_rose_beyond_its_declared_tolerance(summary, baseline) -> None:
    j = (baseline.join(summary.select("level", "rule", pl.col("findings").alias("now")),
                       on=["level", "rule"], how="left")
         .with_columns(pl.col("now").fill_null(0))
         .with_columns((pl.col("now") - pl.col("findings")).alias("delta")))
    over = j.filter(pl.col("delta") > pl.col("tolerance"))
    assert over.height == 0, (
        "a QC finding grew beyond its tolerance. Say in the commit message which chunk brought it and "
        f"why it is acceptable, then re-record rules/qc_advisories.tsv:\n{over}")


def test_a_rule_that_stopped_firing_is_noticed_rather_than_silently_kept(summary, baseline) -> None:
    """A count falling to zero is good news, and it still has to be recorded.

    Otherwise the baseline keeps a stale allowance: `alpha and beta cdr3 identical` was 98 rows until
    #583 cleared them, and a baseline that still allowed 98 would let them come back unnoticed.
    """
    now = {(r["level"], r["rule"]): r["findings"] for r in summary.iter_rows(named=True)}
    stale = [(r["level"], r["rule"], r["findings"]) for r in baseline.iter_rows(named=True)
             if now.get((r["level"], r["rule"]), 0) == 0 and r["findings"] > 0]
    assert not stale, (
        f"these no longer fire and their allowance is stale - remove them from {BASELINE}: {stale}")


def test_a_new_rule_must_be_declared_before_it_can_fire(summary, baseline) -> None:
    known = {(r["level"], r["rule"]) for r in baseline.iter_rows(named=True)}
    new = [(r["level"], r["rule"], r["findings"]) for r in summary.iter_rows(named=True)
           if (r["level"], r["rule"]) not in known]
    assert not new, (
        f"undeclared QC findings: {new}. A new rule firing is either a data defect to fix or a "
        f"baseline row to add with its measured count.")


def test_nothing_fatal_is_hiding_in_the_baseline(baseline) -> None:
    """Only advisories belong here. A fatal rule is gated by `qc --strict` exiting non-zero, and
    giving it a tolerance would turn a hard failure into an allowance."""
    fatal = baseline.filter(~pl.col("advisory"))
    assert fatal.height == 0, f"fatal rules must not carry a tolerance:\n{fatal}"


def _chunk(**cols):
    """A minimal frame with the columns the counter rule reads."""
    import polars as pl

    n = len(next(iter(cols.values())))
    base = {"chunk.file": ["c.txt"] * n, "chunk.row": list(range(1, n + 1)),
            "antigen.epitope": ["SIINFEKL"] * n, "antigen.gene": [""] * n,
            "antigen.species": [""] * n, "mhc.a": [""] * n, "mhc.b": [""] * n}
    return pl.DataFrame(base | cols)


def _counters(frame):
    from vdjdb.qc.rules import COUNTER_COLUMNS, _no_counter

    return {c for c in COUNTER_COLUMNS
            if not frame.select(_no_counter(c).all()).item()}


def test_a_gene_column_holding_a_consecutive_run_is_a_counter():
    """#694: `Eef2`..`Eef188`, one label per row in lockstep with the row index."""
    frame = _chunk(**{"antigen.gene": [f"Eef{n}" for n in range(2, 12)]})
    assert _counters(frame) == {"antigen.gene"}


def test_an_allele_column_holding_a_consecutive_run_is_a_counter():
    """#625: `HLA-A*24:03`..`HLA-A*24:20`, where the paper types every donor `A*24:02`."""
    frame = _chunk(**{"mhc.a": [f"HLA-A*24:{n:02d}" for n in range(3, 21)]})
    assert _counters(frame) == {"mhc.a"}


def test_two_neighbouring_alleles_are_not_a_counter():
    """A paper reporting two restrictions must not trip this - that is #597, a different question."""
    frame = _chunk(**{"mhc.a": ["HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:06"]})
    assert _counters(frame) == set()


def test_a_sparse_family_is_not_a_counter():
    """`TRBV5-1, TRBV5-3, TRBV5-8` fills 3 of an 8-wide span. A counter is dense by construction."""
    frame = _chunk(**{"antigen.gene": ["MAGEA1", "MAGEA3", "MAGEA10"]})
    assert _counters(frame) == set()


def test_one_epitope_does_not_borrow_another_epitope_s_values():
    """The run has to be dense inside one `(chunk, epitope)`, which is where it should be constant."""
    import polars as pl

    frame = pl.DataFrame({
        "chunk.file": ["c.txt"] * 6, "chunk.row": list(range(1, 7)),
        "antigen.epitope": ["AAA", "AAA", "AAA", "BBB", "BBB", "BBB"],
        "antigen.gene": ["Eef2", "Eef2", "Eef2", "Eef3", "Eef3", "Eef3"],
        "antigen.species": [""] * 6, "mhc.a": [""] * 6, "mhc.b": [""] * 6,
    })
    assert _counters(frame) == set(), "constant within each epitope, so nothing is a counter"


def test_the_corpus_carries_exactly_the_one_declared_counter():
    """`PMID_39286976.tsv`, 38 rows, and `patches/mhc.dict` already repairs it at build time.

    The rule starts as a regression guard on a clean corpus, which is the state #597 argues a new
    rule should start from: a finding here means a counter that arrived since, not one of a list a
    reader has to remember to ignore.
    """
    import polars as pl

    baseline = pl.read_csv("rules/qc_advisories.tsv", separator="\t")
    counters = baseline.filter(pl.col("rule").str.starts_with("counter in"))
    assert counters.height == 1
    assert counters["rule"].item() == "counter in mhc.a"
    assert counters["findings"].item() == 38
    assert counters["chunks"].item() == 1
