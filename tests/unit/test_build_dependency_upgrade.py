"""Regression cases for the shared TRA/TRD annotation migration in #1405."""
import polars as pl

from vdjdb.annotate.cdr3fix import markup
from vdjdb.assemble.evidence import structures
from vdjdb.score.evidence import integrate


def test_shared_v_proposals_and_supported_anchor_repairs():
    # One batch across species; these sequences include the conserved junction anchors.
    keys = pl.DataFrame({
        'species': ['HomoSapiens'] * 4 + ['MusMusculus'],
        'cdr3': ['CAMREAQSGGYQKVTF', 'CAFMKDAGGTSYGKLTF', 'SAVSAPGAGANDLFF',
                 'YAENHALHF', 'SAMDSNYQL'],
        'v': ['', 'TRAV38-2/DV8*01', 'TRAV14/DV4*01',
              'TRAV5*01,TRAV23/DV6*01', 'TRAV14-2'],
        'j': ['TRAJ13*01', 'TRAJ52*01', 'TRAJ36*01', 'TRAJ41*01', 'TRAJ33'],
    })
    from vdjdb.curate.nomenclature import harmonise_segments

    harmonised, _ = harmonise_segments(keys.rename({'v': 'v.alpha', 'j': 'j.alpha'}))
    keys = harmonised.rename({'v.alpha': 'v', 'j.alpha': 'j'})
    got = markup(keys, 'alpha')
    expected = [
        ('CAMREAQSGGYQKVTF', 'TRAV14/DV4*01'),
        ('CAFMKDAGGTSYGKLTF', 'TRAV38-1*01'),
        ('CAVSAPGAGANDLFF', 'TRAV8-4*01'),
        ('CAENHALHF', 'TRAV13-2*01'),
        ('CAMDSNYQLIW', 'TRAV13-4/DV7*03'),
    ]
    for source, (junction, v) in zip(keys.iter_rows(named=True), expected, strict=True):
        row = got.filter(pl.col('cdr3') == source['cdr3']).row(0, named=True)
        assert row['v'] == source['v']
        assert (row['__cdr3'], row['__v']) == (junction, v)
        assert row['__gv'] == (v if not source['v'] else '')


def test_recalled_v_associates_structure_and_integrates_existing_evidence():
    records = pl.DataFrame({
        'record_id': ['old', 'recalled'], 'species': ['HomoSapiens'] * 2,
        'antigen.epitope': ['GILGFVFTL'] * 2, 'mhc.a': ['HLA-A*02:01'] * 2,
        'mhc.b': ['B2M'] * 2, 'mhc.class': ['MHCI'] * 2,
        'reference.id': ['paper1', 'paper2'], 'meta.subject.id': ['donor1', 'donor2'],
        'vdjdb.score': [1, 0],
    })
    chains = pl.DataFrame({
        'record_id': ['old', 'old', 'recalled', 'recalled'],
        'gene': ['TRA', 'TRB'] * 2,
        'cdr3': ['CAFMKDAGGTSYGKLTF', 'CASSIGSYGYTF'] * 2,
        'v.segm': ['TRAV38-1*01', 'TRBV19*01', 'TRAV38-2/DV8*01', 'TRBV19*01'],
        'j.segm': ['TRAJ52*01', 'TRBJ1-2*01'] * 2,
        'TCR_hash': ['current1'] * 2 + ['current2'] * 2,
    })
    metadata = records.head(1).drop('record_id').with_columns(
        pl.lit('source').alias('tcr_pmhc_hash'), pl.lit('false').alias('is_native'),
        pl.lit('CAFMKDAGGTSYGKLTF').alias('cdr3.alpha'),
        pl.lit('CASSIGSYGYTF').alias('cdr3.beta'), pl.lit('0.88').alias('iptm'))
    before = structures(chains, metadata, records=records)
    assert before.is_empty(), 'conflicting current V genes must prevent association'
    corrected = chains.with_columns(pl.when(pl.col('gene') == 'TRA')
                                   .then(pl.lit('TRAV38-1*01'))
                                   .otherwise(pl.col('v.segm')).alias('v.segm'))
    evidence = structures(corrected, metadata, records=records)
    assert set(evidence['record_id']) == {'old', 'recalled'}
    integrated = integrate(records, corrected, evidence)
    assert integrated.height == 1
    assert integrated['legacy_score'][0] == 1
    assert integrated['donors'][0] == integrated['a_studies'][0] == 2
    assert integrated['iptm'][0] == 0.88


def test_rhesus_d_posterior_allows_only_machine_rounding():
    import pytest
    from vdjtools.model.reference import translate

    from vdjdb.annotate.junction import infer

    keys = pl.DataFrame({'species': ['MacacaMulatta'], 'cdr3': ['CASSVIRTSNQPQYF'],
                         'v.segm': ['TRBV14*01'], 'j.segm': ['TRBJ1-5*01']})
    row = infer(keys).row(0, named=True)
    assert row['d.inferred'] == 'TRBD1'
    assert row['d.posterior'] == pytest.approx(0.8072042852066905, rel=0, abs=2e-15)
    assert translate(row['cdr3nt']) == keys['cdr3'][0]


def test_import_preserves_caller_openmp_capacity():
    import os
    import subprocess
    import sys

    env = dict(os.environ)
    env.pop('OMP_THREAD_LIMIT', None)
    env.pop('POLARS_MAX_THREADS', None)
    probe = ("import os, vdjdb, vdjtools; "
             "assert 'OMP_THREAD_LIMIT' not in os.environ; "
             "assert 'POLARS_MAX_THREADS' not in os.environ")
    result = subprocess.run([sys.executable, '-c', probe], env=env,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
