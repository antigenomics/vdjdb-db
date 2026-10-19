"""Evidence integration preserves observational and missingness boundaries."""
import polars as pl
import pytest

from vdjdb.score.evidence import categorical, integrate


def test_categories_keep_unknown_separate_from_zero():
    frame = pl.DataFrame({'a_studies': [None, 0, 1, 2, 9], 'iptm': [None, .2, .5, .8, 1.]})
    got = categorical(frame)
    assert got['a_studies'].to_list() == ['?', '0', '1', '2', '3+']
    assert got['iptm'].to_list() == ['?', 'low', 'mid', 'high', 'high']
    with pytest.raises(ValueError):
        categorical(pl.DataFrame({'iptm': [1.2]}))


def test_integrates_built_score_counts_and_motif_missingness():
    records = pl.DataFrame({
        'record_id': ['r1', 'r2', 'r3'], 'species': ['HomoSapiens'] * 3,
        'antigen.epitope': ['PEP'] * 3, 'mhc.a': ['HLA-A*02:01'] * 3,
        'mhc.b': ['B2M'] * 3, 'mhc.class': ['MHCI'] * 3,
        'reference.id': ['P1', 'P2', 'P1'], 'meta.subject.id': ['D1', 'D1', 'D1'],
        'meta.replica.id': ['E1', 'E2', 'E1'], 'meta.study.id': ['S'] * 3,
        'method.identification': ['tetramer-sort,cultured-T-cells',
                                  'phage display, magnetic beads', 'tetramer-sort'],
        'method.verification': ['', 'direct', ''], 'method.singlecell': ['yes'] * 3,
        'method.sequencing': ['sanger'] * 3, 'vdjdb.score': [1, 3, 1],
    })
    chains = pl.DataFrame({
        'record_id': ['r1', 'r1', 'r2', 'r2', 'r3', 'r3'],
        'gene': ['TRA', 'TRB'] * 3, 'cdr3': ['CAVF', 'CASSF', 'CAVF', 'CATSF', 'CAVF', 'CASSF'],
        'v.segm': ['TRAV1*01', 'TRBV1*01', 'TRAV1*02', 'TRBV1*01', 'TRAV1*01', 'TRBV1*01'],
        'j.segm': ['TRAJ1*01', 'TRBJ1*01'] * 3,
    })
    evidence = pl.DataFrame({
        'record_id': ['r1', 'r1', 'r2', 'r1', 'r2'], 'gene': ['TRA', 'TRB', '', 'TRA', 'TRA'],
        'evidence_type': ['motif_tcrnet', 'motif_tcremp_nonmember', 'structure_model',
                          'motif_tcremp_nonmember', 'motif_tcremp'],
        'evidence_value': ['', '', 'iptm', '', ''], 'evidence_score': [5., 0., .9, 0., 4.],
        'evidence_source': ['', '', 'model1', '', ''],
    })
    got = integrate(records, chains, evidence)
    r1 = got.filter(pl.col('b_junction') == 'CASSF').row(0, named=True)
    assert got.height == 2
    assert r1['a_studies'] == 2 and r1['legacy_score'] == 1
    assert r1['donors'] == 1
    assert r1['a_tcrnet'] == 1 and r1['b_tcremp'] == 0
    assert r1['b_tcrnet'] is None and r1['iptm'] is None


def test_nan_is_unknown_and_fractional_counts_are_rejected():
    assert categorical(pl.DataFrame({'iptm': [float('nan')]}))['iptm'][0] == '?'
    with pytest.raises(ValueError):
        categorical(pl.DataFrame({'a_studies': [1.5]}))


def test_training_selected_iptm_bounds_preserve_missing_and_validate_order():
    frame = pl.DataFrame({'iptm': [None, .6, .8, .85, .9]})
    assert categorical(frame, iptm_bounds=(.8, .9))['iptm'].to_list() == ['?', 'low', 'mid', 'mid', 'high']
    for bounds in [(.9, .8), (.5, .5), (float('nan'), .8), (-.1, .8), (.5, 1.1)]:
        with pytest.raises(ValueError):
            categorical(frame, iptm_bounds=bounds)
