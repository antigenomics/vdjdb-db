"""CDR3 markup against external nucleotide truth, on VDJdb's own records.

`assemble.master.fix_cdr3` documented its two engines on **coverage** - how often each declines - and
coverage cannot tell "answers more" from "answers better". The phase 5 swap decision (ROADMAP
section 28) rests on that measurement alone. These tests supply the missing half, and they do it on
VDJdb records rather than on synthetic sequences.

There is one engine now. The k-mer scanner was deleted with `res/` (#658), so its arm here is the
`vEnd` and `jStart` the **2026-06-03 release shipped** rather than a re-run of vendored code. The two
agree to within a fraction of a point (see
`test_the_shipped_scanner_answers_on_almost_every_truth_row`), and a release cannot drift from what
shipped while a vendored copy can.

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
shipped `v.segm` and `j.segm`, which carry an allele on 99.2 % of cells, and the allele used to change
the answer: arda declined `v.end` for any allele whose `cdr3_anchors.tsv` status is `truncated`, where
the same gene without a suffix resolved to `*01` and mapped. That was antigenomics/arda#135, opened
from this fixture and fixed in arda 2.31.0, which places the boundary a truncated germline supports
and marks it `TruncatedGermline`. Both decline assertions below are now plain assertions, and on this
overlap arda declines nothing.

These assertions are written to hold both before and after that arda release, so the pin can move
without a test rewrite. What they refuse is a *regression*: an engine that declines more than it
used to, or that walks further into the junction than the scanner VDJdb already ships.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.assemble.master import _markup_arda
from vdjdb.convert.coords import nt_to_aa_boundary_expr
from vdjdb.emit.vdjdb3 import read_table

pytestmark = pytest.mark.release

CONTROL = Path.home() / "hf/airr_control/human.trb.ntvj.vdjtools.tsv.gz"

#: The released zip, read for the **scanner's own answer**. The k-mer scanner was deleted with
#: `res/` (#658), so the comparison arm is no longer a re-run of vendored code: it is the `vEnd` and
#: `jStart` the 2026-06-03 release actually shipped, which is what that code produced and is frozen.
#: Better than re-running it, because a vendored copy can drift from what shipped and a release
#: cannot.
REFERENCE = Path(os.environ.get("VDJDB_REFERENCE_ZIP", "reference.zip"))

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
def truth_all() -> pl.DataFrame:
    """VDJdb human TRB records whose junction the control also observed, with an agreed boundary."""
    directory = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    needed = [directory / f"{n}.parquet" for n in ("records", "chains")]
    if missing := [str(p) for p in needed if not p.exists()]:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    if not CONTROL.exists():
        pytest.skip(f"{CONTROL} is not checked out; see CLAUDE.md developer setup")
    records, chains = (read_table(directory, n) for n in ("records", "chains"))
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
    got = _markup_arda(only, "beta").select(
        *KEY, pl.col("__vend").alias("arda.v_end"), pl.col("__jstart").alias("arda.j_start"))
    out = out.join(got, on=KEY, how="left", maintain_order="left")
    return out.join(_shipped_scanner(), left_on=["cdr3", "vg", "jg"],
                    right_on=["cdr3", "vg", "jg"], how="left", maintain_order="left")


@pytest.fixture(scope="module")
def truth(truth_all) -> pl.DataFrame:
    """Compare engines only where the released database contains the sequence/gene key."""
    return truth_all.filter(pl.col("legacy.present").fill_null(False))


@pytest.mark.parametrize("coord", ["v_end", "j_start"])
def test_current_engine_accuracy_includes_new_imports(truth_all, coord) -> None:
    near = truth_all.filter(
        (pl.col(f"arda.{coord}") >= 0)
        & ((pl.col(f"arda.{coord}") - pl.col(f"t.{coord}")).abs() <= 1)).height
    assert near / truth_all.height > 0.95


def _shipped_scanner() -> pl.DataFrame:
    """The scanner's `vEnd` / `jStart` as the 2026-06-03 release shipped them, per gene-level key.

    Keyed the same way the control is - `(cdr3, V gene, J gene)` - because the reference's allele
    suffix is the one the scanner resolved and the nomenclature phase has since corrected 10,605 of
    them, so an allele-level key would lose exactly the rows this comparison is about. 123,165
    distinct beta keys, of which 45 carry two `vEnd` values and 386 two `jStart`; their boundaries are null
    rather than averaged; their keys remain present, on the same rule the control's own boundary uses.
    """
    if not REFERENCE.exists():
        pytest.skip(f"{REFERENCE} is not present; set VDJDB_REFERENCE_ZIP")
    from vdjdb.compare.diff import Bundle, _read_table

    ref = _read_table(Bundle(REFERENCE).read_bytes("vdjdb.txt"))
    return (ref.filter(pl.col("gene") == "TRB")
            .select("cdr3",
                    pl.col("v.segm").str.split("*").list.first().alias("vg"),
                    pl.col("j.segm").str.split("*").list.first().alias("jg"),
                    pl.col("cdr3fix").str.json_path_match("$.vEnd").cast(pl.Int64).alias("v_end"),
                    pl.col("cdr3fix").str.json_path_match("$.jStart").cast(pl.Int64)
                      .alias("j_start"))
            .group_by("cdr3", "vg", "jg")
            .agg(pl.col("v_end").n_unique().alias("nv"), pl.col("j_start").n_unique().alias("nj"),
                 pl.col("v_end").first().alias("legacy.v_end"),
                 pl.col("j_start").first().alias("legacy.j_start"))
            .with_columns(
                pl.when((pl.col("nv") == 1) & (pl.col("nj") == 1))
                  .then(pl.col("legacy.v_end")).otherwise(None).alias("legacy.v_end"),
                pl.when((pl.col("nv") == 1) & (pl.col("nj") == 1))
                  .then(pl.col("legacy.j_start")).otherwise(None).alias("legacy.j_start"),
                pl.lit(True).alias("legacy.present"))
            .drop("nv", "nj"))


def test_the_overlap_is_large_enough_to_draw_a_conclusion_from(truth) -> None:
    """Stated as a test so a shrinking reference is a failure, not a quietly weaker number."""
    assert truth.height >= MIN_RECORDS
    assert truth["cdr3"].n_unique() >= MIN_RECORDS // 2


def test_the_shipped_scanner_answers_on_almost_every_truth_row(truth) -> None:
    """The comparison arm is the release's own cells now, so a **null** means the release has no row
    for that key - not that the scanner declined. Every head-to-head below filters on `>= 0`, which
    drops a null silently, so the null count is what has to be gated: a join that quietly stopped
    matching would read as "the scanner agrees with everything".

    34 of 1,632 on the 2026-06-03 release. Measured against a re-run of the vendored scanner before it
    was deleted, the two arms agree closely - `v.end` exact 72.90 % here against 72.86 % re-run,
    `j.start` 95.81 % against 97.06 % - which is what makes reading the release faithful.
    """
    missing = truth.filter(pl.col("legacy.v_end").is_null()).height
    assert missing <= truth.height // 10, (
        f"the release has no boundary for {missing} of {truth.height} truth rows; the gene-level "
        f"join has stopped matching, so every head-to-head below is measuring a smaller set than it "
        f"reports")


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
    """Why the swap was proposed, now holding on both coordinates.

    On `j.start` arda declines strictly less often, which is the half that was never in doubt and is
    the larger half: the swap took `j.start` coverage from 277,939 to 284,880 of 286,047 chains.

    On `v.end` it used to decline *more*, and the cause was never the alignment. `cdr3_anchors.tsv`
    marks 63 human V alleles `status = truncated`, because IMGT ships those allele records as partial
    sequences that stop inside the anchor region, and `cdr3fix` answered `FailedBadSegment` for all
    of them - including the 38 whose `templated_aa` is still 3 residues or longer and places a
    boundary perfectly well. That was `antigenomics/arda#135`, opened from this fixture's own
    measurement and fixed in arda 2.31.0, which places the boundary those germlines support and marks
    it `TruncatedGermline` so the caller can see it is a lower bound.

    Measured on the current corpus, 285,989 chains: `v.end` unmapped falls from 5,307 to 4,163, and
    1,144 chains move `FailedBadSegment` -> `TruncatedGermline` with `fix.good` false -> true on
    1,139 of them. `TRBV11-2*02` is 850 of the 1,144. No chain loses a boundary and no `cdr3` changes.
    On this fixture's overlap arda now declines `v.end` on 0 rows, against 0 for the scanner.
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


def test_the_fallback_is_unexercised_here_because_arda_declines_on_nothing(truth) -> None:
    """Why the bound above now checks no rows, asserted rather than left to be discovered.

    The bound is vacuous on an empty set and `continue` is exactly how it would go quiet, so the
    emptiness needs a stated cause. It has one: arda 2.31.0 places a boundary on every row of this
    overlap, both coordinates, so there is nothing left for a fallback to fill here. Before it, 7
    rows exercised the `v.end` fallback (4 exact, 7 within one residue).

    This is the notification if that reverses. arda declining again makes the count non-zero, this
    test fails, and the bound above starts checking rows on the same run.

    The fallback is still filled in the build - 2,442 chains for `v.end.inferred` and 488 for
    `j.start.inferred` - on chains this control does not observe. `tests/release/test_tables_contract.py`
    gates those counts; only the comparison against external nucleotides is what has run out of rows.
    """
    for coord, fallback in (("v_end", "v.end.inferred"), ("j_start", "j.start.inferred")):
        declined = truth.filter(pl.col(f"arda.{coord}") < 0)
        assert declined.height == 0, (
            f"arda declines {coord} on {declined.height} truth rows, so the fallback bound above is "
            f"checking them again - read its result rather than this test's")
        assert truth.filter(pl.col(fallback) >= 0).height == 0, (
            f"{fallback} is filled on a truth row where arda answered; the fallback is supposed to "
            f"be masked to -1 wherever the markup engine mapped the boundary")
