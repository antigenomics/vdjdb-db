"""Junction-nucleotide inference (#461): call resolution, determinism, and back-translation."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.annotate import junction


class FakeModel:
    """Just the ``genomic`` tables :func:`junction._resolver` reads."""

    def __init__(self) -> None:
        self.genomic = {
            "genes_v": pl.DataFrame({
                "gene": ["TRBV1", "TRBV2", "TRBV2"],
                "v_allele": ["TRBV1*01", "TRBV2*01", "TRBV2*02"],
            }),
            "genes_j": pl.DataFrame({"gene": ["TRBJ1"], "j_allele": ["TRBJ1*01"]}),
        }


def test_a_known_allele_is_passed_through():
    r = junction._resolver(FakeModel(), "genes_v", "v_allele")
    assert r["TRBV2*01"] == "TRBV2*01"


def test_a_gene_with_one_modelled_allele_resolves_to_it():
    """No invention: there is exactly one thing `TRBV1` can mean to this model."""
    assert junction._resolver(FakeModel(), "genes_v", "v_allele")["TRBV1"] == "TRBV1*01"


def test_a_gene_with_several_alleles_marginalises_rather_than_guessing():
    """Picking one would be the #327 mistake -- an allele asserted from nothing."""
    assert junction._resolver(FakeModel(), "genes_v", "v_allele").get("TRBV2") is None


def test_an_unknown_or_empty_call_marginalises():
    r = junction._resolver(FakeModel(), "genes_v", "v_allele")
    assert r[""] is None
    assert r.get("TRBV 2-7") is None      # a space where a dash belongs: real, and phase 9's
    assert r.get("TRBV99*01") is None


# -- against the real models -------------------------------------------------------------------

KEYS = pl.DataFrame({
    "cdr3": ["CASSIRSSYEQYF", "CASSLAPGATNEKLFF", "CASSPGQGAYEQYF", "CASSQDRGNTGELFF"],
    "v.segm": ["TRBV10-3*01", "TRBV7-9*01", "TRBV5-1*01", "TRBV4-1*01"],
    "j.segm": ["TRBJ2-7*01", "TRBJ1-4*01", "TRBJ2-7*01", "TRBJ2-2*01"],
})


def test_the_inferred_nucleotides_back_translate_to_the_junction_they_came_from():
    """#461's acceptance criterion. 0 mismatches over the whole corpus."""
    from vdjtools.model import translate

    got = junction.infer(KEYS, "HomoSapiens", "TRB", workers=1)
    resolved = got.filter(pl.col("cdr3nt") != "")
    assert resolved.height == KEYS.height
    for aa, nt in zip(resolved["cdr3"], resolved["cdr3nt"], strict=True):
        assert translate(nt) == aa


def test_the_worker_count_does_not_change_the_answer():
    """CLAUDE.md rule 7. Contiguous slices over a sorted key set, reassembled in slice order."""
    one = junction.infer(KEYS, "HomoSapiens", "TRB", workers=1)
    four = junction.infer(KEYS, "HomoSapiens", "TRB", workers=4)
    assert one.equals(four)


def test_a_species_with_no_model_gets_empty_columns_rather_than_an_error():
    """The corpus holds 1,402 MacacaMulatta chains and no macaque model exists."""
    got = junction.infer(KEYS, "MacacaMulatta", "TRB")
    assert got["cdr3nt"].to_list() == [""] * KEYS.height
    assert got["cdr3nt.pgen"].null_count() == KEYS.height


@pytest.mark.parametrize("missing", ["cdr3", "v.segm", "j.segm"])
def test_a_chain_missing_a_call_keeps_its_row_and_gets_no_nucleotides(missing):
    chains = pl.DataFrame({
        "record_id": ["r1", "r2"], "gene": ["TRB", "TRB"],
        "cdr3": ["CASSIRSSYEQYF", "CASSPGQGAYEQYF"],
        "v.segm": ["TRBV10-3*01", "TRBV5-1*01"], "j.segm": ["TRBJ2-7*01", "TRBJ2-7*01"],
    }).with_columns(pl.when(pl.col("record_id") == "r2").then(pl.lit(""))
                    .otherwise(pl.col(missing)).alias(missing))
    records = pl.DataFrame({"record_id": ["r1", "r2"],
                            "species": ["HomoSapiens", "HomoSapiens"]})
    got = junction.add_junction_nt(chains, records, workers=1)
    assert got.height == 2, "a chain is never dropped for being un-inferrable"
    assert got.filter(pl.col("record_id") == "r2")["cdr3nt"][0] == ""
    assert got.filter(pl.col("record_id") == "r1")["cdr3nt"][0] != ""


def test_the_worker_count_actually_buys_wall_time():
    """CLAUDE.md section 0e: a pool that never ran is indistinguishable from a slow one.

    The correctness test above passes whether or not the threads run in parallel, so it cannot see a
    pool that silently serialised. This measures it. `infer_nt` calls into vdjtools' native path,
    which releases the GIL for part of the work, so the speedup is real and sublinear. Measured on a
    16-core M3, four workers: **2.26x over these 400 keys**, 2.04x over 1,200 and 2.20x over 3,000, so
    the key count is not what sets it.

    **The bar depends on the core count, because the quantity does.** Four worker threads on a 4-vCPU
    GitHub runner share those four cores with polars' own pool and the interpreter, and there is no
    headroom to win: the same 400 keys measured **1.29x** there on 2026-09-28 and failed a flat 1.3x
    bar that had only ever been checked on a 16-core laptop. A bar that fails on a correct pool is
    worse than a loose one, because the next person deletes it.

    So: 1.7x where there are at least 8 cores, which is a much tighter detector than the 1.3x it
    replaces and still 25 % below the measured 2.26x; 1.15x on 4 to 7, which is what a 4-vCPU runner
    can show while still being nowhere near the 1.0x a dead pool gives. The parallel leg is timed
    twice and the faster run counts, because a shared runner's scheduler noise is one-sided -- a live
    pool can be unlucky, a dead one never gets faster on a retry.
    """
    import itertools
    import os
    import time

    from vdjdb.timing import cores_available

    # The cores this process may use, not the ones the machine has. `os.cpu_count()` reports 40 on a
    # SLURM `-c 4` cpuset, which broke this test in two directions at once on 2026-09-28: it chose the
    # `>= 8` bar of 1.7x for an allocation that had four cores, and it compared a load average of 33.9
    # against 40 rather than 4, so the guard below that exists for exactly that case did not fire.
    # Four workers bought 1.56x and the test failed on a pool that was working correctly.
    cores = cores_available()
    if cores < 4:
        pytest.skip("needs at least 4 cores to say anything about scaling")
    # ⚠ Refuse to measure a machine that cannot be measured. A ratio taken while the box is
    # oversubscribed is not a property of the pool: measured 2.26x quiet, 1.62x inside a full suite
    # run and 1.48x at load average 22 on 16 cores, all with the same correct pool. Skipping is the
    # only outcome the measurement supports - the alternative is a bar loose enough to pass under any
    # load, which is a bar that no longer detects a dead pool (ROADMAP_local section 57.4).
    #
    # The load average is the whole machine's even under a cpuset, which is the behaviour wanted here:
    # a four-core slice of a node that forty other processes are using cannot be timed either.
    load = os.getloadavg()[0]
    if load > cores:
        pytest.skip(f"load average {load:.1f} on {cores} cores: too busy to time anything")

    aa = "ACDEFGHIKLMNPQRSTVWY"
    # Distinct junctions with the V and J anchors of a real key untouched, so every one resolves: the
    # timing would mean nothing if half of them took the no-scenario fast path. Asserted below.
    mids = ["".join(p) for p in itertools.product(aa, repeat=2)]
    keys = pl.DataFrame({"cdr3": [f"CASS{m}IRSSYEQYF" for m in mids],
                         "v.segm": ["TRBV10-3*01"] * len(mids),
                         "j.segm": ["TRBJ2-7*01"] * len(mids)}).sort("cdr3", "v.segm", "j.segm")

    junction.infer(keys.head(20), "HomoSapiens", "TRB", workers=1)   # load the model, once

    def timed(workers: int) -> tuple[float, object]:
        start = time.perf_counter()
        out = junction.infer(keys, "HomoSapiens", "TRB", workers=workers)
        return time.perf_counter() - start, out

    # Best of two on **both** legs. One-sided noise on either leg moves the ratio, and the machine
    # this runs on is not quiet: measured 2.26x in isolation and 1.62x inside a full suite run at load
    # average 19 on 16 cores, which failed a 1.7x bar on a pool that was working correctly. Taking the
    # best of two on each side compares best case with best case, which is what a ratio under load
    # means (ROADMAP_local section 57.4).
    serial, serial_out = min((timed(1) for _ in range(2)), key=lambda r: r[0])
    parallel = min(timed(4)[0] for _ in range(2))
    bar = 1.7 if cores >= 8 else 1.15

    assert (serial_out["cdr3nt"] != "").all(), "the timing set must not take the no-scenario path"
    assert serial / parallel > bar, (
        f"four workers bought {serial / parallel:.2f}x over {keys.height} keys "
        f"({serial:.2f} s -> {parallel:.2f} s) on {cores} cores at load {load:.1f}, "
        f"bar {bar}x; "
        f"the pool is overhead rather than parallelism")
