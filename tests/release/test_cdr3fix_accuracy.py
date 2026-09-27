"""Is `arda` more *accurate* than the legacy k-mer scanner, or only more willing to answer?

`assemble.master.fix_cdr3` documents the two engines on **coverage** - how often each declines -
and coverage cannot tell "answers more" from "answers better". The phase 5 swap decision (ROADMAP
section 28) rests on that measurement alone, so these tests supply the missing half.

Neither engine is scored against the other. Both are scored against a reference that knows the
answer for a reason independent of either alignment strategy:

* `model` - `vdjtools.model.infer_nt` finds the most probable recombination scenario for a junction
  from OLGA-learned statistics. `annotate.junction._infer_slice` already computes it for every chain
  in the build and keeps only the D coordinates, so `v_end` / `j_start` are free here. Not an oracle:
  a third method, but an independent one.
* `nucleotides` - `airr_control`'s `human.trb.ntvj` carries real clonotypes with the **observed**
  `cdr3nt` and `VEnd` / `JStart` in nucleotide space, derived where the germline boundary is an
  alignment fact rather than a guess from a protein. Neither engine sees the nucleotides.
  CC-BY-NC-ND, read as an input, aggregate statistics only, nothing derived from it is written
  (hard rule 5). Skipped when the checkout is absent.

Measured 2026-09-27, both references agreeing on direction: the legacy is **more concordant**, arda
is **more complete**, and a legacy-first / arda-fallback hybrid beats both. The assertions below are
directional with slack rather than frozen counts, so an arda release that changes the alignment
shows up without a fixture rewrite every time a chunk lands.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.assemble.master import _markup_arda, _markup_legacy

pytestmark = pytest.mark.release

SEED = 20260927

#: Small enough that both engines run in seconds, large enough that a several-hundred-row
#: head-to-head gap cannot be sampling noise.
N = 4_000

#: `v.end` and `j.start` are documented as "0-based amino acid, junction space" and their
#: closedness is not stated. Fitted in both directions against both references: each is the
#: **count of residues the segment touches**, so a nucleotide index converts by ceiling, not by
#: `coords.nt_to_aa`, which floors. Getting this wrong moves every rate by tens of percent and
#: changes no head-to-head count, which is why the assertions below are head-to-head.
def _to_aa(nt: pl.Expr) -> pl.Expr:
    return (nt + 2) // 3


CONTROL = Path.home() / "hf/airr_control/human.trb.ntvj.vdjtools.tsv.gz"


KEY = ["species", "cdr3", "v", "j"]


def _mark(keys: pl.DataFrame) -> pl.DataFrame:
    """Mark up with both engines and join the result beside whatever else `keys` carries.

    `_markup_legacy` unpacks `iter_rows()` into exactly four names, so it is handed a frame of
    exactly the four key columns; the truth columns come back on the join.
    """
    out = keys
    only = keys.select(KEY)
    for name, run in (("legacy", _markup_legacy), ("arda", _markup_arda)):
        got = run(only, "beta").select(
            *KEY, pl.col("__vend").alias(f"{name}.v_end"),
            pl.col("__jstart").alias(f"{name}.j_start"))
        out = out.join(got, on=KEY, how="left", maintain_order="left")
    return out


def _head_to_head(m: pl.DataFrame, coord: str, truth: str) -> tuple[int, int, int]:
    """(legacy alone right, arda alone right, rows both engines map)."""
    both = m.filter((pl.col(f"legacy.{coord}") >= 0) & (pl.col(f"arda.{coord}") >= 0))
    only_l = both.filter((pl.col(f"legacy.{coord}") == pl.col(truth))
                         & (pl.col(f"arda.{coord}") != pl.col(truth))).height
    only_a = both.filter((pl.col(f"arda.{coord}") == pl.col(truth))
                         & (pl.col(f"legacy.{coord}") != pl.col(truth))).height
    return only_l, only_a, both.height


# -- reference 1: the recombination model, on VDJdb's own records ------------------------------

@pytest.fixture(scope="module")
def against_model() -> pl.DataFrame:
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    needed = [directory / f"{n}.parquet" for n in ("records", "chains")]
    if missing := [str(p) for p in needed if not p.exists()]:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    records, chains = (pl.read_parquet(p) for p in needed)
    keys = (chains.join(records.select("record_id", "species"), on="record_id", how="inner")
            .filter((pl.col("species") == "HomoSapiens") & (pl.col("gene") == "TRB")
                    & (pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))
            .select("species", "cdr3", pl.col("v.segm").alias("v"), pl.col("j.segm").alias("j"))
            .unique(maintain_order=True).sort("species", "cdr3", "v", "j"))
    sample = keys.sample(min(N, keys.height), seed=SEED).sort("species", "cdr3", "v", "j")

    from vdjtools.model import infer_nt, load_bundled

    from vdjdb.annotate.junction import _resolver
    model = load_bundled("TRB", "olga", organism="human")
    vmap = _resolver(model, "genes_v", "v_allele")
    jmap = _resolver(model, "genes_j", "j_allele")
    ve, js = [], []
    for cdr3, v, j in sample.select("cdr3", "v", "j").iter_rows():
        s = infer_nt(model, cdr3, v=vmap.get(v), j=jmap.get(j))
        ve.append(None if s is None else s.v_end)
        js.append(None if s is None else s.j_start)
    ref = (sample.with_columns(pl.Series("nt.v_end", ve, dtype=pl.Int64),
                               pl.Series("nt.j_start", js, dtype=pl.Int64))
           .filter(pl.col("nt.v_end").is_not_null()))
    # `Scenario.v_end` is half-open where the control's `VEnd` is closed (measured: +1 nt on 65 %
    # of real clonotypes), so the closed form is what converts the same way as the control.
    return _mark(ref.select("species", "cdr3", "v", "j",
                            _to_aa(pl.col("nt.v_end") - 1).alias("t.v_end"),
                            _to_aa(pl.col("nt.j_start")).alias("t.j_start")))


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_the_legacy_agrees_with_the_recombination_model_more_often_than_arda(against_model, coord):
    """Both engines are given the same junction and the same V/J calls; one third method judges.

    The gap is one-sided, not a wash: on the rows both engines answer, the legacy is right where
    arda is wrong far more often than the reverse.
    """
    only_l, only_a, both = _head_to_head(against_model, coord, f"t.{coord}")
    assert both > N // 2, "too few rows mapped by both engines to compare"
    assert only_l > only_a, (
        f"{coord}: legacy alone right {only_l}, arda alone right {only_a} on {both} rows - "
        "if this has reversed, the phase 5 swap decision should be revisited with the new numbers")
    assert only_l >= 3 * max(only_a, 1), (
        f"{coord}: the legacy's lead has narrowed to {only_l}:{only_a}; it was measured at "
        "382:1 (v_end) and 912:0 (j_start) on 2026-09-27")


def test_arda_declines_far_less_often_than_the_legacy(against_model):
    """The other half, and the reason the swap was ever proposed: arda almost always answers."""
    counts = {e: against_model.filter(pl.col(f"{e}.v_end") < 0).height for e in ("legacy", "arda")}
    assert counts["arda"] < counts["legacy"], counts
    assert counts["arda"] <= 0.2 * counts["legacy"], (
        f"arda used to decline a small fraction of what the legacy declines: {counts}")


# -- reference 2: nucleotide-level markup of real repertoires ----------------------------------

@pytest.fixture(scope="module")
def against_nucleotides() -> pl.DataFrame:
    if not CONTROL.exists():
        pytest.skip(f"{CONTROL} is not checked out; see CLAUDE.md developer setup")
    d = pl.read_csv(CONTROL, separator="\t", n_rows=400_000,
                    schema_overrides={"cdr3nt": pl.Utf8, "cdr3aa": pl.Utf8, "v": pl.Utf8,
                                      "j": pl.Utf8})
    d = (d.filter(pl.col("v").str.starts_with("TRBV") & pl.col("j").str.starts_with("TRBJ")
                  & ~pl.col("v").str.contains(",") & ~pl.col("j").str.contains(",")
                  & (pl.col("VEnd") >= 0) & (pl.col("JStart") >= 0)
                  & pl.col("cdr3aa").str.starts_with("C") & pl.col("cdr3aa").str.contains("[FW]$")
                  & (pl.col("cdr3nt").str.len_chars() == 3 * pl.col("cdr3aa").str.len_chars()))
         .unique(subset=["cdr3aa", "v", "j"], keep="first", maintain_order=True)
         .sample(N, seed=SEED))
    keys = d.select(pl.lit("HomoSapiens").alias("species"), pl.col("cdr3aa").alias("cdr3"),
                    (pl.col("v") + "*01").alias("v"), (pl.col("j") + "*01").alias("j"),
                    _to_aa(pl.col("VEnd")).alias("t.v_end"),
                    _to_aa(pl.col("JStart")).alias("t.j_start"))
    return _mark(keys)


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_real_nucleotide_markup_agrees_with_the_legacy_more_often_too(against_nucleotides, coord):
    """The independent reference that is not a model: the nucleotides were actually observed."""
    only_l, only_a, both = _head_to_head(against_nucleotides, coord, f"t.{coord}")
    assert both > N // 2
    assert only_l > only_a, f"{coord}: legacy {only_l} vs arda {only_a} on {both} rows"


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_legacy_first_with_arda_as_fallback_beats_either_engine_alone(against_nucleotides, coord):
    """What the two halves imply: take the more concordant answer, fall back for coverage.

    Neither engine dominates - the legacy is right more often, arda answers more often - and the
    combination is better than both on the same rows, with nothing left unmapped. This is the
    measurement the swap decision needs, so it is pinned rather than left in a working note.
    """
    m, t = against_nucleotides, f"t.{coord}"
    hybrid = pl.when(pl.col(f"legacy.{coord}") >= 0).then(pl.col(f"legacy.{coord}")) \
               .otherwise(pl.col(f"arda.{coord}"))
    got = m.with_columns(hybrid.alias("hybrid"))
    exact = {name: got.filter(pl.col(name) == pl.col(t)).height
             for name in (f"legacy.{coord}", f"arda.{coord}", "hybrid")}
    assert exact["hybrid"] >= exact[f"legacy.{coord}"], exact
    assert exact["hybrid"] >= exact[f"arda.{coord}"], exact
    assert got.filter(pl.col("hybrid") < 0).height == 0, "the fallback must leave nothing unmapped"
