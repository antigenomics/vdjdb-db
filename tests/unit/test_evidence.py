"""Independent-study evidence: what counts as replication, and what does not."""
from __future__ import annotations

import polars as pl

from vdjdb.assemble.evidence import build_evidence, independent_study, support_counts
from vdjdb.schema import EVIDENCE_TABLE_COLUMNS


def tables(*rows: tuple[str, str, str, str]) -> tuple[pl.DataFrame, pl.DataFrame]:
    """``(record_id, clonotype, epitope, reference)`` -> the two frames evidence reads."""
    records = pl.DataFrame({
        "record_id": [r[0] for r in rows],
        "antigen.epitope": [r[2] for r in rows],
        "reference.id": [r[3] for r in rows],
    })
    chains = pl.DataFrame({
        "record_id": [r[0] for r in rows],
        "gene": ["TRB"] * len(rows),
        "clonotype_id": [r[1] for r in rows],
    })
    return records, chains


def test_one_paper_reporting_a_clonotype_twice_is_not_replication():
    """Two records, one reference: the same study, so nothing is independent about it."""
    ev = independent_study(*tables(("r1", "c", "GIL", "PMID:1"), ("r2", "c", "GIL", "PMID:1")))
    assert ev.is_empty()


def test_two_papers_reporting_one_clonotype_epitope_pair_is_evidence():
    records, chains = tables(("r1", "c", "GIL", "PMID:1"), ("r2", "c", "GIL", "PMID:2"))
    ev = independent_study(records, chains)
    assert ev.height == 2
    assert ev["evidence_score"].to_list() == [2.0, 2.0]
    # the evidence a reader wants is who *else* saw it, not the paper they are already reading
    assert dict(zip(ev["record_id"], ev["evidence_value"], strict=True)) == {
        "r1": "PMID:2", "r2": "PMID:1"}


def test_the_same_clonotype_against_a_different_epitope_is_not_support():
    """Support is per clonotype-epitope pair: cross-reactivity is not replication."""
    ev = independent_study(*tables(("r1", "c", "GIL", "PMID:1"), ("r2", "c", "NLV", "PMID:2")))
    assert ev.is_empty()


def test_a_different_clonotype_in_the_same_epitope_is_not_support():
    ev = independent_study(*tables(("r1", "c1", "GIL", "PMID:1"), ("r2", "c2", "GIL", "PMID:2")))
    assert ev.is_empty()


def test_support_counts_and_the_evidence_rows_are_the_same_computation():
    """ROADMAP §11.1 tunes on this; the shipped column must not be a second implementation."""
    records, chains = tables(("r1", "c", "GIL", "PMID:1"), ("r2", "c", "GIL", "PMID:2"),
                             ("r3", "c", "GIL", "PMID:2"), ("r4", "d", "GIL", "PMID:9"))
    counts = support_counts(records, chains)
    supported = set(counts.filter(pl.col("studies") > 1)["clonotype_id"])
    assert supported == {"c"}
    ev = independent_study(records, chains)
    assert set(ev["record_id"]) == {"r1", "r2", "r3"}
    assert ev["evidence_score"].unique().to_list() == [2.0]


def test_the_evidence_key_is_unique_and_the_column_order_is_the_declared_one():
    records, chains = tables(("r1", "c", "GIL", "PMID:1"), ("r2", "c", "GIL", "PMID:2"))
    ev = build_evidence(records, chains, release="v1")
    assert tuple(ev.columns) == EVIDENCE_TABLE_COLUMNS
    assert ev.select("record_id", "evidence_id").n_unique() == ev.height
    assert ev["first_seen_release"].unique().to_list() == ["v1"]


def test_no_producer_yields_an_empty_table_with_the_right_schema():
    """A build with nothing to say must still write a readable evidence table."""
    ev = build_evidence(*tables(("r1", "c", "GIL", "PMID:1")))
    assert ev.is_empty()
    assert tuple(ev.columns) == EVIDENCE_TABLE_COLUMNS
