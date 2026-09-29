"""IMGT's F / ORF / P verdict on a named segment (#634).

`imgt_alleles.tsv.gz` is a committed, reviewed authority and `functionality` was read by nothing:
`vdjdb qc` asks whether a call *looks* like a TRBV name and `curate.nomenclature` asks whether IMGT
*has* it, and neither asks whether IMGT thinks the gene is functional.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate import functionality as fn


@pytest.mark.parametrize(("verdict", "want"), [
    ("F", True), ("(F)", True), ("[F]", True),
    ("ORF", False), ("(ORF)", False), ("P", False), ("[P]", False),
    # A joined verdict is a *gene*'s, and the question a finding answers is whether the call is
    # evidence of a problem - so one functional allele is enough.
    ("F/ORF", True), ("(F)/F", True), ("ORF/P", False),
    ("", False), ("nonsense", False),
])
def test_a_verdict_is_functional_when_any_allele_is(verdict: str, want: bool) -> None:
    """`any`, not `all`, and the difference was 15,283 findings.

    Under `all`, `TRBJ2-7` reads `F/ORF` because `*02` is an ORF, and it is one of the commonest J
    calls in VDJdb - the advisory rule fired on 17,891 chunk rows instead of 2,608 chain-segments, and
    nothing a curator can act on was in the difference.
    """
    assert fn.is_functional(verdict) is want


def test_the_verdict_is_species_specific() -> None:
    """IMGT has `TRBV7-1*01` as ORF in human and P in *Macaca fascicularis*.

    #634 read that as one allele carrying two verdicts and asked for a tie-break rule. It is entirely
    cross-species: keyed on `(species, allele)` the table has **zero** duplicate rows, so the
    resolution is a lookup and not a reconciliation.
    """
    table = fn.verdicts()
    assert table[("HomoSapiens", "TRBV7-1*01")] == ("ORF", "allele")
    assert len({k for k in table}) == len(list(table)), "one verdict per (species, call)"


def test_an_allele_verdict_beats_its_gene() -> None:
    """Allele first: an allele-level verdict answers the question asked, a gene-level one infers it."""
    table = fn.verdicts()
    assert table[("HomoSapiens", "TRAJ58*01")][1] == "allele"
    assert table[("HomoSapiens", "TRAJ24")][1] == "gene"


def _chains(*calls: tuple[str, str]) -> tuple[pl.DataFrame, pl.DataFrame]:
    chains = pl.DataFrame(
        [{"record_id": f"r{i}", "gene": "TRB", "v.segm": v, "j.segm": j}
         for i, (v, j) in enumerate(calls)])
    records = pl.DataFrame({"record_id": [f"r{i}" for i in range(len(calls))],
                            "species": ["HomoSapiens"] * len(calls)})
    return chains, records


def test_only_a_non_functional_call_is_reported() -> None:
    flagged = fn.report(*_chains(("TRBV10-3*01", "TRBJ2-7*01"),      # F, F
                                ("TRBV21-1*01", "TRBJ2-7*01"),      # P on the V
                                ("TRBV10-3*01", "TRAJ58*01")))      # ORF on the J
    assert set(flagged["call"]) == {"TRBV21-1*01", "TRAJ58*01"}
    assert flagged.filter(call="TRBV21-1*01")["functionality"][0] == "P"
    assert tuple(flagged.columns) == fn.COLUMNS


def test_a_call_imgt_lists_nowhere_is_not_this_check_s_finding() -> None:
    """`curate.nomenclature` is the check for a name IMGT does not have. Two checks, two findings."""
    flagged = fn.report(*_chains(("TRBV999*01", "TRBJ2-7*01")))
    assert flagged.is_empty()


def test_a_blank_call_is_not_a_finding() -> None:
    assert fn.report(*_chains(("", "TRBJ2-7*01"))).is_empty()


def test_the_summary_counts_chains_per_verdict() -> None:
    flagged = fn.report(*_chains(("TRBV21-1*01", "TRBJ2-7*01"),
                                ("TRBV21-1*01", "TRBJ2-7*01"),
                                ("TRBV10-3*01", "TRAJ58*01")))
    got = fn.summarise(flagged)
    assert got.filter(segment="V", functionality="P")["chains"][0] == 2
    assert got["chains"].sum() == flagged.height


def test_an_empty_report_summarises_to_an_empty_frame_with_the_right_shape() -> None:
    """A frame whose schema depends on whether it has rows is the pandas ambiguity returning."""
    got = fn.summarise(fn.report(*_chains(("TRBV10-3*01", "TRBJ2-7*01"))))
    assert got.is_empty()
    assert "chains" in got.columns
