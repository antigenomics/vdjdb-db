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


def test_support_references_are_individual_values():
    records, chains = tables(('r1','c','GIL','PMID:1'),('r2','c','GIL','PMID:2'),
                             ('r3','c','GIL','PMID:3'))
    ev = independent_study(records, chains)
    assert ev.height == 6
    assert ev.filter(pl.col('record_id')=='r1')['evidence_value'].to_list() == ['PMID:2','PMID:3']
    assert ev.select('record_id','evidence_id').n_unique() == 6
    assert ev['evidence_score'].unique().to_list() == [3.]


def test_structures_link_exact_complexes_and_keep_metrics_separate():
    from vdjdb.assemble.evidence import structures
    chains = pl.DataFrame({'record_id':['r1','r1','r2','r3'], 'gene':['TRA','TRB','TRA','TRA'],
                           'TCR_hash':['h1','h1','h1','h2']})
    metadata = pl.DataFrame({'tcr_pmhc_hash':['h1','h1','h1','unmatched'],
        'is_native':['False','False','True','False'], 'meta.structure.id':['','','1abc',''],
        'iptm':['0.7','0.7','1.0','0.9'], 'num_contacts':['10','10','12','15']})
    ev = structures(chains,metadata,release='v1')
    assert ev['record_id'].unique().sort().to_list() == ['r1','r2']
    assert ev['gene'].unique().to_list() == ['']
    assert ev.height == 6  # native + two model metrics, once per record, never per chain
    scores = ev.filter(pl.col('evidence_type')=='structure_model')['evidence_score'].sort().to_list()
    assert scores == [.7,.7,10.,10.]
    assert ev.select('record_id','evidence_id').n_unique() == ev.height
    import pytest
    with pytest.raises(ValueError,match='conflicting'):
        structures(chains,pl.concat([metadata, metadata.head(1).with_columns(pl.lit('0.8').alias('iptm'))]))


def test_motif_evidence_distinguishes_nonmembership_from_unavailable():
    from vdjdb.assemble.evidence import motifs
    records = pl.DataFrame({'record_id':['r1','r2','r3','r4'], 'species':['HomoSapiens']*4,
                            'antigen.epitope':['P','P','Q','P']})
    chains = pl.DataFrame({'record_id':['r1','r2','r3','r4'],'gene':['TRA']*4,
        'cdr3':['CASF','CATF','CASF','CASF'],'v.segm':['V1','V1','V1','V2'],'j.segm':['J1']*4})
    cohort = pl.DataFrame({'species':['HomoSapiens']*2,'gene':['TRA']*2,
        'antigen.epitope':['P']*2,'junction_aa':['CASF','CATF'],'v_call':['V1']*2,'j_call':['J1']*2})
    members = cohort.head(1).with_columns(pl.lit('cluster1').alias('cid'),pl.lit(5).alias('csz'))
    ev = motifs(records,chains,cohort,members,method='tcrnet',release='v1')
    assert dict(zip(ev['record_id'],ev['evidence_score'], strict=True)) == {'r1':5.,'r2':0.,'r4':5.}
    assert ev.filter(pl.col('record_id')=='r2')['evidence_type'].item() == 'motif_tcrnet_nonmember'
    emp = motifs(records,chains,cohort,members,method='tcremp')
    assert emp['record_id'].to_list() == ['r1','r2']
