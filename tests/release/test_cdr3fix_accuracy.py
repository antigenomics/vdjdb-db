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

The `declines` column of that table was measured on gene-level V/J calls. This fixture passes the
shipped `v.segm` and `j.segm`, which carry an allele on 99.2 % of cells, and the allele changes the
answer: arda declines `v.end` for any allele whose `cdr3_anchors.tsv` status is `truncated`, where
the same gene without a suffix resolves to `*01` and maps. That is antigenomics/arda#135, and it is
why the `v.end` decline assertion below is a strict xfail while the `j.start` one is not.

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
from vdjdb.convert.coords import nt_to_aa_boundary_expr

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
#:
#: It was spelled out here and nowhere else until the model boundary started shipping as a fallback
#: (#631). `vdjdb.convert.coords` owns it now, in both a scalar and a column form that
#: `tests/unit/test_coords.py` asserts agree - two spellings of one coordinate conversion is what
#: `CLAUDE.md`'s four-spaces table warns about.
_aa = nt_to_aa_boundary_expr


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
                    "v.end", "j.start", "v.end.inferred", "j.start.inferred",
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


@pytest.mark.xfail(strict=True, reason="antigenomics/arda#135: arda declines v_end on alleles "
                                       "whose cdr3_anchors.tsv status is truncated")
def test_arda_declines_far_less_often_than_the_shipped_scanner(truth) -> None:
    """Why the swap was proposed - and the one place it does not hold, which is a V-side allele gap.

    On `j.start` arda declines strictly less often, which is the half that was never in doubt and is
    the larger half: the swap took `j.start` coverage from 277,939 to 284,880 of 286,047 chains.

    On `v.end` it declines *more*, and the cause is not the alignment. `cdr3_anchors.tsv` marks 63
    human V alleles `status = truncated`, because IMGT ships those allele records as partial
    sequences that stop inside the anchor region, and `cdr3fix` answers `FailedBadSegment` for all
    of them - including the 38 whose `templated_aa` is still 3 residues or longer and would place a
    boundary perfectly well. Measured 2026-09-28 over the 114,117 distinct human TRB
    `(cdr3, v.segm, j.segm)` keys in the built corpus: the scanner declines `v.end` on 1,791, arda
    on 2,672, and 2,566 keys are declined by arda while the scanner maps them. 1,702 of those carry
    no allele at all and mostly name a family rather than a gene (`TRBV6`, `TRBV12`), where
    declining is the correct answer; 851 carry an allele, 837 of them `*02`, and `TRBV11-2*02` alone
    is 802.

    Kept as a strict xfail rather than relaxed, because the 7 cases this fixture can check against
    external nucleotides say the boundary is genuinely there: truth places `v.end` at 4 or 5, the
    scanner matches it exactly on 5 of 7 and within one residue on 7 of 7, and arda answers -1 on
    all 7. So there is nothing to concede here - the assertion is right and arda#135 is the bug.
    Delete the marker when that ships; strict makes the pass itself the notification.
    """
    counts = {e: truth.filter(pl.col(f"{e}.v_end") < 0).height for e in ("legacy", "arda")}
    assert counts["arda"] <= counts["legacy"], counts


def test_arda_declines_j_start_less_often_than_the_shipped_scanner(truth) -> None:
    """The half that does hold, asserted separately so arda#135 cannot mask a J-side regression."""
    counts = {e: truth.filter(pl.col(f"{e}.j_start") < 0).height for e in ("legacy", "arda")}
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


# -- the model boundary as a fallback where arda declines (#631) --------------------------------

def test_the_model_boundary_is_never_more_than_one_residue_out_where_arda_declines(truth) -> None:
    """The fallback's accuracy against the external truth, on the rows it is actually used on.

    ⚠ **The comparable set is 7 rows, not the 55 #631 states.** That figure does not reproduce at this
    overlap: `truth` needs a human TRB record whose junction the 400,000-row control also observed
    *and* on which every control observation agrees about the boundary, and only 7 of those are rows
    where arda declined and the model answered. Measured on all 7: **4 exact (0.571) and 7 within one
    residue (1.000)**.

    So the rate is recorded and the *bound* is gated. An exact-match bar on n = 7 is a bar on noise,
    and a bar that fails on correct code is one the next person deletes (`ROADMAP_local.md` §52.2).
    Off-by-at-most-one is a property rather than a rate: it says the model is picking the right codon
    neighbourhood, which is the claim that justifies shipping it where `-1` carries nothing at all.
    """
    for coord, fallback, engine in (("v_end", "v.end.inferred", "arda.v_end"),
                                    ("j_start", "j.start.inferred", "arda.j_start")):
        used = truth.filter((pl.col(engine) < 0) & (pl.col(fallback) >= 0))
        if not used.height:
            continue          # no row in the overlap needs the fallback for this coordinate
        off = (used[fallback] - used[f"t.{coord}"]).abs()
        assert int(off.max()) <= 1, (
            f"{coord}: the fallback is {int(off.max())} residues out on "
            f"{int((off > 1).sum())} of {used.height} rows where arda declined\n"
            f"{used.filter(off > 1).select('cdr3', 'v', 'j', fallback, f't.{coord}').head(5)}")


def test_the_comparable_set_for_the_fallback_has_not_vanished(truth) -> None:
    """Stated as a test so a shrinking overlap is a failure rather than a silently weaker number.

    The bound above is vacuous on an empty set, and `continue` is exactly how it would go quiet.
    """
    v_used = truth.filter((pl.col("arda.v_end") < 0) & (pl.col("v.end.inferred") >= 0))
    assert v_used.height >= 5, (
        f"only {v_used.height} truth rows exercise the v.end fallback; it was 7 on 2026-09-29, so "
        f"either the overlap shrank or arda stopped declining")
