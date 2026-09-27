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
    which releases the GIL for part of the work, so the speedup is real and sublinear: measured on a
    16-core M3, 1.52x at two workers, **2.27x at the shipped four**, 2.87x at eight, and the same
    2.28x on 3,000 real corpus keys rather than these 400 synthetic ones.

    The bar is 1.3x, well below the measured 2.27x and well above the 1.0x a dead pool would give.
    Loose on purpose: this is a regression detector for the pool disappearing, not a benchmark.
    """
    import itertools
    import os
    import time

    if (os.cpu_count() or 1) < 4:
        pytest.skip("needs at least 4 cores to say anything about scaling")

    aa = "ACDEFGHIKLMNPQRSTVWY"
    # Distinct junctions with the V and J anchors of a real key untouched, so every one resolves: the
    # timing would mean nothing if half of them took the no-scenario fast path. Asserted below.
    mids = ["".join(p) for p in itertools.product(aa, repeat=2)]
    keys = pl.DataFrame({"cdr3": [f"CASS{m}IRSSYEQYF" for m in mids],
                         "v.segm": ["TRBV10-3*01"] * len(mids),
                         "j.segm": ["TRBJ2-7*01"] * len(mids)}).sort("cdr3", "v.segm", "j.segm")

    junction.infer(keys.head(20), "HomoSapiens", "TRB", workers=1)   # load the model, once
    serial_start = time.perf_counter()
    serial_out = junction.infer(keys, "HomoSapiens", "TRB", workers=1)
    serial = time.perf_counter() - serial_start
    parallel_start = time.perf_counter()
    junction.infer(keys, "HomoSapiens", "TRB", workers=4)
    parallel = time.perf_counter() - parallel_start

    assert (serial_out["cdr3nt"] != "").all(), "the timing set must not take the no-scenario path"
    assert serial / parallel > 1.3, (
        f"four workers bought {serial / parallel:.2f}x over {keys.height} keys "
        f"({serial:.2f} s -> {parallel:.2f} s); the pool is overhead rather than parallelism")
