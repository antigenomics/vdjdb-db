"""The retired build's two master-table checks, replayed over the whole corpus and accounted for.

`tests/unit/test_legacy_qc_parity.py` proves the two builds agree on one constructed record per
failure mode. This proves it on all 192,793 records, which is the half a fixture cannot do: a mode
nobody thought to construct still has to be accounted for here.

``runBuidDatabase.py`` ran three checks after CDR3 repair and filed what failed into
``vdjdb_full_gene_broken.txt``, ``vdjdb_full_allele_broken.txt`` and ``vdjdb_full_cdr3aa_broken.txt``.
Nothing in the rewrite replaced the first two: ``build_master`` calls ``harmonise_segments`` and
discards its report, and that report names what was *rewritten* rather than what could not be, so a
call naming a gene no authority carries reached every shipped table with no report anywhere. The third
was replaced by something better and, until the ``unanchored`` fallback landed, also narrower.

So each of the three is **partitioned with no remainder**: every row it flagged is either flagged by
the report that replaced it, or falls in a class this file names and checks against the reference. A
count comparison would pass while the two sets drifted apart; a partition cannot.

Measured 2026-09-29 over the built corpus:

=====================  ============  ================================================  ============
legacy finding         chains/calls  accounted for by                                  chains/calls
=====================  ============  ================================================  ============
`cdr3 not C..[WF]`     956 chains    `anchors.tsv`                                     475
                                     the J germline really does end in that residue    481
`gene not in IMGT`     2,582 calls   `nomenclature.tsv`                                2,582
`allele out of range`  1 call        IMGT lists the allele                             1
=====================  ============  ================================================  ============

The 481 are legacy false positives and the reason the hardcoded rule had to go: ``TRAJ35*01``
templates ``IGFGNVLHC`` and mouse ``TRAJ7*01`` templates ``DYSNNRLTL``, so "ends in Phe or Trp" calls
269 and 212 correct junctions broken. The one allele finding is the same kind of mistake in the other
direction: ``alleles_match_check`` compared ``int(allele)`` against a per-gene count from a 741-row
human immunoglobulin table, so it was a range check rather than a membership one, and IMGT does list
``TRBV28*02``.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

import legacy_qc as L
from vdjdb.curate.anchors import templated
from vdjdb.curate.nomenclature import _SPLIT, _imgt
from vdjdb.emit.vdjdb3 import to_internal

pytestmark = pytest.mark.release

#: Legacy's `is_qq_seq_biologically_valid` on the shipped junctions. A rise means a chunk landed with
#: junctions in the wrong coordinate space; the partition below is what says whether that is reported.
LEGACY_JUNCTION_FINDINGS = 956
#: Legacy's `gene_match_check`, human only, as the driver ran it.
LEGACY_GENE_FINDINGS = 2582
#: A band, because both numbers move with every chunk. What must not move is the *remainder*, which is
#: asserted at zero.
TOLERANCE = 400


def _root() -> Path:
    return Path(os.environ.get("VDJDB_ROOT", "."))


@pytest.fixture(scope="module")
def chains() -> pl.DataFrame:
    """The shipped chains with their species, which is what both legacy checks read."""
    tables = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    missing = [n for n in ("chains.tsv", "records.tsv") if not (tables / n).exists()]
    if missing:
        pytest.skip(f"no built tables: {', '.join(missing)}")
    # The TSV, not the parquet: these checks assert on what shipped as text. `to_internal` puts the
    # shipped `underscore_case` header back under the names the rest of this file uses.
    ch = to_internal(pl.read_csv(tables / "chains.tsv", separator="\t", infer_schema_length=0))
    rec = to_internal(pl.read_csv(tables / "records.tsv", separator="\t", infer_schema_length=0))
    return ch.join(rec.select("record_id", "species"), on="record_id", how="left")


@pytest.fixture(scope="module")
def reports() -> dict[str, pl.DataFrame]:
    directory = Path(os.environ.get("VDJDB_REPORTS", "out/reports"))
    missing = [n for n in ("anchors.tsv", "nomenclature.tsv") if not (directory / n).exists()]
    if missing:
        pytest.skip(f"no build reports: {', '.join(missing)}")
    return {n: pl.read_csv(directory / n, separator="\t", infer_schema_length=0)
            for n in ("anchors.tsv", "nomenclature.tsv")}


# --------------------------------------------------------------------------------------------
# vdjdb_full_cdr3aa_broken.txt -> anchors.tsv
# --------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def retired_junction_findings(chains) -> pl.DataFrame:
    """Every chain the retired `is_qq_seq_biologically_valid` rejects.

    A Python call per row, on purpose. The predicate is simple enough to write as a polars expression
    and writing it as one would make this a reimplementation of the oracle rather than a use of it,
    which is the only thing that makes the partition below evidence. 0.3 s over 285,989 chains, and
    module-scoped so it is paid once.
    """
    return chains.filter(
        ~pl.col("cdr3").map_elements(L.is_qq_seq_biologically_valid, return_dtype=pl.Boolean))


def _reported(anchors: pl.DataFrame) -> set[tuple[str, str]]:
    return set(zip(anchors["record_id"], anchors["gene"], strict=True))


def test_every_junction_the_retired_check_rejected_is_reported_or_justified(
        retired_junction_findings, reports) -> None:
    """The partition with no remainder. A junction legacy called broken is either reported by the
    germline-based check, or the germline it names really does end in that residue."""
    flagged = retired_junction_findings
    assert flagged.height == pytest.approx(LEGACY_JUNCTION_FINDINGS, abs=TOLERANCE), (
        f"the retired junction check now rejects {flagged.height} chains against a recorded "
        f"{LEGACY_JUNCTION_FINDINGS}. Re-measure and say in the commit which chunk moved it.")

    reported = _reported(reports["anchors.tsv"])
    remainder = flagged.filter(
        ~pl.struct("record_id", "gene").map_elements(
            lambda s: (s["record_id"], s["gene"]) in reported, return_dtype=pl.Boolean))

    # The germline that justifies each one, read out of the reference rather than assumed.
    justified = remainder.with_columns(
        pl.struct("species", "j.segm").map_elements(
            lambda s: (templated(s["species"] or "", "J", s["j.segm"] or "") or "")[-1:],
            return_dtype=pl.Utf8).alias("anchor"))
    unaccounted = justified.filter(
        (pl.col("anchor") == "") | (pl.col("cdr3").str.slice(-1) != pl.col("anchor")))
    assert unaccounted.height == 0, (
        f"{unaccounted.height} junction(s) the retired build rejected are reported by nothing and "
        "justified by no germline:\n"
        f"{unaccounted.select('record_id', 'gene', 'species', 'cdr3', 'j.segm', 'anchor').head(10)}")


# A global reported fraction depends on the mix of valid germline-supported junctions.
# The per-record partition above is the gate: no unexplained rejection may disappear.


def test_a_junction_with_no_germline_to_check_is_still_reported(reports) -> None:
    """The `unanchored` fallback, which is the whole reason the partition closes.

    With no J call, or a call the reference does not have, there is no germline to compare against and
    the check used to decline silently - 118 chains that the retired build's crude rule was the only
    thing covering. The universal anchor is the fallback and the finding says so.
    """
    unanchored = reports["anchors.tsv"].filter(pl.col("defect").str.contains("unanchored"))
    assert unanchored.height > 0, (
        "no chain is reported as `unanchored`. Either the corpus has none, which would be new, or the "
        "fallback in `anchors.classify` has stopped firing and 118 chains went silent.")
    # No germline means nothing to compare against, so no sibling allele can be named either -
    # naming one would be inventing a call. (The *sequence* repair is `arda.cdr3fix`'s and this
    # module no longer computes one at all, #711.)
    assert unanchored.filter(pl.col("sibling.call") != "").height == 0, (
        "a sibling allele was named for a junction with no germline behind it")


# --------------------------------------------------------------------------------------------
# vdjdb_full_gene_broken.txt -> nomenclature.tsv
# --------------------------------------------------------------------------------------------

def test_every_segment_call_the_retired_check_rejected_is_reported(chains, reports) -> None:
    """Partitioned the same way: reported, or a multi-call whose every part IMGT has.

    The second class is legacy's own defect. A curator recording two candidates writes
    ``TRBD1,TRBD2``, which is not an allele name, and ``gene_match_check`` split on ``*`` alone and so
    rejected the whole string.
    """
    root = _root()
    tables = _imgt(root)
    reported = set(zip(reports["nomenclature.tsv"]["species"], reports["nomenclature.tsv"]["call"],
                       strict=True))

    flagged, unaccounted = 0, []
    for column in ("v.segm", "j.segm"):
        for species, call, n in chains.group_by("species", column).len().rows():
            if not call or species != "HomoSapiens" or L.gene_match_check(call, root):
                continue
            flagged += n
            if (species, call) in reported:
                continue
            alleles, genes = tables[species]
            parts = [p for p in _SPLIT.split(call.strip()) if p]
            if all(p in alleles or p.split("*")[0] in genes for p in parts):
                continue        # a multi-call legacy could not read; every part is a real IMGT name
            unaccounted.append((species, column, call, n))

    assert not unaccounted, (
        f"segment calls the retired build rejected that no report names: {unaccounted}")
    assert flagged == pytest.approx(LEGACY_GENE_FINDINGS, abs=TOLERANCE), (
        f"the retired gene check now rejects {flagged} chain-calls against a recorded "
        f"{LEGACY_GENE_FINDINGS}. Re-measure `nomenclature.tsv` and say which chunk moved it.")


def test_the_one_allele_finding_the_retired_check_produced_is_a_false_positive(chains) -> None:
    """`alleles_match_check` was a range check wearing the clothes of a membership check.

    It compared ``int(allele)`` against a per-gene count, so ``*07`` passed and ``*08`` failed whether
    or not IMGT listed either. Its single finding on the whole corpus is an allele IMGT has, which is
    why nothing replaced it: `nomenclature.tsv` asks the membership question directly.
    """
    root = _root()
    alleles, _ = _imgt(root)["HomoSapiens"]
    human = chains.filter(pl.col("species") == "HomoSapiens")
    for column in ("v.segm", "j.segm"):
        for call in human[column].unique().drop_nulls():
            if not call or L.alleles_match_check(call, root):
                continue
            assert call in alleles, (
                f"{call} fails the retired allele range check and IMGT does not list it either, so "
                "it is a real finding and `nomenclature.tsv` has to carry it")


def test_the_nomenclature_report_separates_the_two_reasons_a_call_lands_in_it(reports) -> None:
    """A family name and a name with no IMGT candidate need different curation, so the report has to
    tell them apart. ``TRBV6`` is nine genes the record did not choose between; ``TRBV28-0`` is not a
    name at all."""
    r = reports["nomenclature.tsv"]
    family = r.filter(pl.col("family.members").cast(pl.Int64) > 0)
    assert family.height > 0 and family.height < r.height, (
        "the report no longer separates an under-specified family from an unrecognised name")
    assert (family["candidates"] != "").all(), "a family finding names no candidates"
