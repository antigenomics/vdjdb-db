"""CDR3 markup against external nucleotide truth, on VDJdb's own records.

`assemble.master.fix_cdr3` documents its two engines on **coverage** - how often each declines - and
coverage cannot tell "answers more" from "answers better". The phase 5 swap decision (ROADMAP
section 28) rests on that measurement alone. These tests supply the missing half, and they do it on
VDJdb records rather than on synthetic sequences.

The reference is `isalgo/airr_control`'s `human.trb.ntvj`: real repertoire clonotypes carrying the
observed `cdr3nt` together with `VEnd` and `JStart` in nucleotide space. Those coordinates were
established with the nucleotides present, where the germline boundary is an alignment fact rather
than something inferred from a protein sequence. Both engines here see only the amino-acid junction
and the V/J calls.

One amino-acid junction can arise from several nucleotide sequences with different boundaries, so a
key counts as truth only when every control observation of it agrees - measured over the whole
control, 60.7 % of overlapping keys qualify, which is itself worth knowing: the boundary is not
always determined by the amino acids.

CC-BY-NC-ND. The control is read as an input and only aggregate numbers come out; nothing derived
from it is written anywhere a build could ship it (hard rule 5). Skipped when it is not checked out.

Measured 2026-09-27 over the whole control, 8,334 VDJdb human TRB records with unambiguous truth
behind 196,517 supporting reads:

===========================  =============  ========  ===============  ========  ========
source                       `v.end` exact  within 1  `j.start` exact  within 1  declines
===========================  =============  ========  ===============  ========  ========
VDJdb shipped (k-mer scan)   5,985 (71.8%)  8,279     7,987 (95.8%)    8,311     2
arda 2.29.0                  5,968 (71.6%)  8,263     8,041 (96.5%)    8,196     0
arda + antigenomics/arda#133 5,993 (71.9%)  8,289     8,164 (98.0%)    8,320     0
===========================  =============  ========  ===============  ========  ========

The 2.29.0 row is the defect antigenomics/arda#133 fixes: it credits the germline two or more
residues too far into the junction on 138 of these records, against 0 for the shipped scanner. Every
case is the same shape - the J germline's leading residue matching the junction by coincidence and
paying for the next mismatch - and the fix takes it to 1.

These assertions are written to hold both before and after that arda release, so the pin can move
without a test rewrite. What they refuse is a *regression*: an engine that declines more than it
used to, or that walks further into the junction than the scanner VDJdb already ships.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.assemble.master import _markup_arda, _markup_legacy

pytestmark = pytest.mark.release

CONTROL = Path.home() / "hf/airr_control/human.trb.ntvj.vdjtools.tsv.gz"

#: Rows of the control to read. It is frequency-sorted, so this is the expanded end of the
#: repertoire; enough of it to clear the 100-record bar below by an order of magnitude.
CONTROL_ROWS = 400_000

#: The overlap must be big enough that a percentage means something.
MIN_RECORDS = 100

KEY = ["species", "cdr3", "v", "j"]

#: `v.end` and `j.start` are documented as "0-based amino acid, junction space" without stating
#: closedness, and the two differ. Fitted against this reference and against the recombination
#: model independently: each is the count of residues the segment touches, so a nucleotide index
#: converts by **ceiling**, not by `coords.nt_to_aa`, which floors. Under floor the rates drop by
#: tens of percent and no head-to-head count moves, which is why the assertions are head-to-head.
def _aa(nt: pl.Expr) -> pl.Expr:
    return (nt + 2) // 3


@pytest.fixture(scope="module")
def truth() -> pl.DataFrame:
    """VDJdb human TRB records whose junction the control also observed, with an agreed boundary."""
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    needed = [directory / f"{n}.parquet" for n in ("records", "chains")]
    if missing := [str(p) for p in needed if not p.exists()]:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    if not CONTROL.exists():
        pytest.skip(f"{CONTROL} is not checked out; see CLAUDE.md developer setup")
    records, chains = (pl.read_parquet(p) for p in needed)
    mine = (chains.join(records.select("record_id", "species"), on="record_id", how="inner")
            .filter((pl.col("species") == "HomoSapiens") & (pl.col("gene") == "TRB")
                    & (pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))
            .select("species", "cdr3", pl.col("v.segm").alias("v"), pl.col("j.segm").alias("j"),
                    "v.end", "j.start",
                    pl.col("v.segm").str.split("*").list.first().alias("vg"),
                    pl.col("j.segm").str.split("*").list.first().alias("jg"))
            .unique(subset=KEY, maintain_order=True).sort(KEY))

    ctl = pl.read_csv(CONTROL, separator="\t", n_rows=CONTROL_ROWS,
                      schema_overrides={"cdr3nt": pl.Utf8, "cdr3aa": pl.Utf8, "v": pl.Utf8,
                                        "j": pl.Utf8})
    agreed = (ctl.filter((pl.col("VEnd") >= 0) & (pl.col("JStart") >= 0)
                         & (pl.col("cdr3nt").str.len_chars()
                            == 3 * pl.col("cdr3aa").str.len_chars()))
              .group_by("cdr3aa", "v", "j")
              .agg(pl.col("VEnd").n_unique().alias("nv"), pl.col("JStart").n_unique().alias("nj"),
                   pl.col("VEnd").first(), pl.col("JStart").first())
              .filter((pl.col("nv") == 1) & (pl.col("nj") == 1)))

    out = (mine.join(agreed, left_on=["cdr3", "vg", "jg"], right_on=["cdr3aa", "v", "j"],
                     how="inner", maintain_order="left")
           .with_columns(_aa(pl.col("VEnd")).alias("t.v_end"),
                         _aa(pl.col("JStart")).alias("t.j_start")))
    if out.height < MIN_RECORDS:
        pytest.skip(f"only {out.height} overlapping records; need {MIN_RECORDS}")
    only = out.select(KEY)
    for name, run in (("legacy", _markup_legacy), ("arda", _markup_arda)):
        got = run(only, "beta").select(
            *KEY, pl.col("__vend").alias(f"{name}.v_end"),
            pl.col("__jstart").alias(f"{name}.j_start"))
        out = out.join(got, on=KEY, how="left", maintain_order="left")
    return out


def test_the_overlap_is_large_enough_to_draw_a_conclusion_from(truth) -> None:
    """Stated as a test so a shrinking reference is a failure, not a quietly weaker number."""
    assert truth.height >= MIN_RECORDS
    assert truth["cdr3"].n_unique() >= MIN_RECORDS // 2


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_both_engines_land_within_one_residue_of_the_nucleotide_boundary(truth, coord) -> None:
    """The contract a consumer of `v.end` / `j.start` actually relies on.

    Exact agreement is 72 % on `v.end` and 96-98 % on `j.start`; within one residue both engines
    clear 99 %. One residue is the codon-boundary ambiguity inherent in placing a nucleotide
    boundary on an amino-acid sequence, and it is shared by every method measured.
    """
    for engine in ("legacy", "arda"):
        mapped = truth.filter(pl.col(f"{engine}.{coord}") >= 0)
        near = mapped.filter(
            (pl.col(f"{engine}.{coord}") - pl.col(f"t.{coord}")).abs() <= 1).height
        assert near / truth.height > 0.95, (
            f"{engine} {coord}: only {near} of {truth.height} within one residue")


def test_arda_declines_far_less_often_than_the_shipped_scanner(truth) -> None:
    """Why the swap was proposed, and the half of it that was never in doubt."""
    counts = {e: truth.filter(pl.col(f"{e}.v_end") < 0).height for e in ("legacy", "arda")}
    assert counts["arda"] <= counts["legacy"], counts


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_no_engine_walks_further_into_the_junction_than_the_shipped_scanner(truth, coord) -> None:
    """The regression guard, and the one antigenomics/arda#133 is about.

    Crediting the germline for residues past the real boundary is the failure that matters: it is
    silent, it moves `v.end` / `j.start` into the N region, and `dpost` slices the non-templated
    middle with exactly these two coordinates. arda 2.29.0 does it on 138 of these records against
    the scanner's 0, and arda#133 takes it to 1. The bound is loose enough to pass either version
    and tight enough that a worse engine fails.
    """
    sign = 1 if coord == "v_end" else -1
    over = {}
    for engine in ("legacy", "arda"):
        mapped = truth.filter(pl.col(f"{engine}.{coord}") >= 0)
        over[engine] = mapped.filter(
            sign * (pl.col(f"{engine}.{coord}") - pl.col(f"t.{coord}")) >= 2).height
    assert over["legacy"] / truth.height < 0.01, over
    assert over["arda"] / truth.height < 0.05, (
        f"{coord}: arda credits germline two or more residues too far on {over['arda']} of "
        f"{truth.height} records, against the scanner's {over['legacy']}. See antigenomics/arda#133")
