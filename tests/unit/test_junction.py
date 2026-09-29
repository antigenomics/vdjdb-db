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

    got = junction.infer(KEYS, "HomoSapiens", "TRB")
    resolved = got.filter(pl.col("cdr3nt") != "")
    assert resolved.height == KEYS.height
    for aa, nt in zip(resolved["cdr3"], resolved["cdr3nt"], strict=True):
        assert translate(nt) == aa


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
        # The markup engine's boundaries. `build_chains` always produces them, and
        # `add_junction_nt` masks its own fallback against them (#631).
        "v.end": [3, 4], "j.start": [6, 7],
    }).with_columns(pl.when(pl.col("record_id") == "r2").then(pl.lit(""))
                    .otherwise(pl.col(missing)).alias(missing))
    records = pl.DataFrame({"record_id": ["r1", "r2"],
                            "species": ["HomoSapiens", "HomoSapiens"]})
    got = junction.add_junction_nt(chains, records)
    assert got.height == 2, "a chain is never dropped for being un-inferrable"
    assert got.filter(pl.col("record_id") == "r2")["cdr3nt"][0] == ""
    assert got.filter(pl.col("record_id") == "r1")["cdr3nt"][0] != ""


def test_the_batch_call_is_the_per_row_loop_and_is_much_faster_than_it():
    """CLAUDE.md rule 3 and section 0e: batch, and ship the check that says the batch happened.

    Two assertions, because each catches a different way of losing this. **Identical** catches a batch
    entry point that is not the same computation - `infer_nt_batch` documents that its result is the
    per-row loop field for field at any thread count, and this is what holds it to that. **Faster**
    catches a regression to the loop, or a batch call that quietly loops internally: neither changes an
    answer, so nothing else here would notice, and this stage was 87 % of the assembly step.

    **The bar is 1.5x, and it is set from the runner rather than from a laptop.** Measured 2026-09-29,
    16 cores: the ratio barely moves with batch size - 8.69x at 400 keys, 8.46x at 1,200, 9.25x at
    3,000 - and moves a great deal with threads, falling to 3.68x / 3.38x / 3.59x at `threads=4`. So
    the shape of the win is thread count, not batch size, and a 4-vCPU runner should see about 3.5x.

    It does not. CI measured **1.58x** on 400 keys at a reported load average of 0.6, which is the
    thing a load check cannot see: a shared cloud vCPU loses time to steal and a cold cache, and
    neither shows up in `getloadavg`. The previous 2x bar was calibrated on this laptop's 10.5x and
    failed a correct build, which is the failure mode `ROADMAP_local` section 52.2 is about.

    So: **best of three**, because contention is one-sided - it only ever makes a run slower, so the
    fastest of several is the better estimator of what the code does - and a bar of 1.5x, which clears
    the worst single-shot observation and is still half the distance to the 1.0x a de-batched call
    gives. A call that loops internally cannot reach it.

    The machine is checked before it is timed, for the reason the process-pool version of this test
    had to learn: a ratio taken while the box is oversubscribed is not a property of the code
    (ROADMAP_local section 57.4).
    """
    import itertools
    import os
    import time

    from vdjtools.model import infer_nt, infer_nt_batch, load_bundled

    from vdjdb.timing import cores_available

    cores = cores_available()
    load = os.getloadavg()[0]
    if load > cores:
        pytest.skip(f"load average {load:.1f} on {cores} cores: too busy to time anything")

    aa = "ACDEFGHIKLMNPQRSTVWY"
    # Distinct junctions with the V and J anchors of a real key untouched, so every one resolves: the
    # timing would mean nothing if half of them took the no-scenario fast path. Asserted below.
    mids = ["".join(p) for p in itertools.product(aa, repeat=2)]
    cdr3 = [f"CASS{m}IRSSYEQYF" for m in mids]
    model = load_bundled("TRB", "olga", organism="human")
    v, j = ["TRBV10-3*01"] * len(cdr3), ["TRBJ2-7*01"] * len(cdr3)

    infer_nt(model, cdr3[0], v=v[0], j=j[0])          # warm the native path, once

    serial_s, batch_s = [], []
    for _ in range(3):
        t = time.perf_counter()
        serial = [infer_nt(model, c, v=vv, j=jj) for c, vv, jj in zip(cdr3, v, j, strict=True)]
        serial_s.append(time.perf_counter() - t)
        t = time.perf_counter()
        batch = infer_nt_batch(model, cdr3, v=v, j=j)
        batch_s.append(time.perf_counter() - t)

    assert all(s is not None for s in serial), "the timing set must not take the no-scenario path"
    assert [s.cdr3_nt for s in serial] == batch["cdr3_nt"].to_list(), (
        "the batched call must be the per-row loop field for field")
    ratio = min(serial_s) / min(batch_s)
    assert ratio > 1.5, (
        f"batching bought {ratio:.2f}x over {len(cdr3)} keys "
        f"({min(serial_s):.2f} s -> {min(batch_s):.2f} s, best of 3) on {cores} cores at load "
        f"{load:.1f}, bar 1.5x; either the loop is back or the batch call is looping internally")


# -- the model boundary is a fallback, never an override (#631) ---------------------------------

def _two_chains(v_end: list[int], j_start: list[int]) -> tuple[pl.DataFrame, pl.DataFrame]:
    chains = pl.DataFrame({
        "record_id": ["r1", "r2"], "gene": ["TRB", "TRB"],
        "cdr3": ["CASSIRSSYEQYF", "CASSLGQAYEQYF"],
        "v.segm": ["TRBV10-3*01", "TRBV7-9*01"], "j.segm": ["TRBJ2-7*01", "TRBJ2-7*01"],
        "v.end": v_end, "j.start": j_start,
    })
    records = pl.DataFrame({"record_id": ["r1", "r2"],
                            "species": ["HomoSapiens", "HomoSapiens"]})
    return chains, records


def test_the_model_boundary_is_dropped_where_the_alignment_answered():
    """The whole safety property of #631: a fallback cannot become an override.

    `v.end` and `j.start` are what `vdjdb-web` reads out of the `cdr3fix` JSON and what the legacy
    tables carry, so changing one of those cells is a data change belonging to a curation decision.
    Masking here rather than at the emitter is what makes that impossible instead of merely unlikely.
    """
    got = junction.add_junction_nt(*_two_chains([3, 4], [6, 7])).sort("record_id")
    assert got["v.end.inferred"].to_list() == [junction.UNMAPPED] * 2
    assert got["j.start.inferred"].to_list() == [junction.UNMAPPED] * 2
    assert got["v.end"].to_list() == [3, 4], "the alignment's answer is untouched"


def test_the_model_boundary_survives_where_the_alignment_declined():
    got = junction.add_junction_nt(*_two_chains([junction.UNMAPPED, 4],
                                                [6, junction.UNMAPPED])).sort("record_id")
    assert got["v.end.inferred"][0] > 0, "arda declined the V boundary; the model has one"
    assert got["v.end.inferred"][1] == junction.UNMAPPED
    assert got["j.start.inferred"][0] == junction.UNMAPPED
    assert got["j.start.inferred"][1] > 0
    # And the boundaries are inside the junction they describe, in residues.
    for row in got.iter_rows(named=True):
        n = len(row["cdr3"])
        for col in ("v.end.inferred", "j.start.inferred"):
            assert row[col] == junction.UNMAPPED or 0 <= row[col] <= n, f"{col} {row[col]} of {n}"


def test_a_species_with_no_model_gets_an_unmapped_boundary_not_a_null():
    """-1 is what this coordinate space already reads as "not mapped" (rule 6's spirit)."""
    got = junction.add_junction_nt(*_two_chains([junction.UNMAPPED] * 2, [junction.UNMAPPED] * 2))
    chains, records = _two_chains([junction.UNMAPPED] * 2, [junction.UNMAPPED] * 2)
    records = records.with_columns(pl.lit("MacacaMulatta").alias("species"))
    got = junction.add_junction_nt(chains, records)
    assert got["v.end.inferred"].to_list() == [junction.UNMAPPED] * 2
    assert got["v.end.inferred"].null_count() == 0
