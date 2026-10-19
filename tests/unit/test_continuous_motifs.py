"""Continuous motif ordering, geometric support and training-only decisions."""
import numpy as np
import polars as pl
import pytest

from vdjdb.score.membership import separated
from vdjdb.validate.membership import fit_pair, predict_pair, ranking_metrics, select_threshold


def test_separation_missingness_and_binary_recovery():
    b = np.array([0., 0., 1., 1., np.nan])
    s = separated(b, np.array([0., 1., 0., 1., .7]))
    assert s[:4] == pytest.approx([0, 1/3, 2/3, 1])
    assert np.array_equal(s[:4] >= .5, b[:4])
    assert np.isnan(s[-1])
    with pytest.raises(ValueError):
        separated(np.array([.5]), np.array([.5]))
    with pytest.raises(ValueError):
        separated(np.array([1]), np.array([1.1]))


def test_tied_thresholds_and_early_recall_metrics():
    pytest.importorskip('sklearn')
    y = np.array([1, 0, 1, 0])
    p = np.array([.8, .8, .2, .1])
    t = select_threshold(y, p)
    assert t == .2
    m = ranking_metrics(y, p)
    assert m['ap'] == pytest.approx(7/12)
    assert m['roc_auc'] == pytest.approx(.625)
    # Tied top group has precision 1/2; no order-dependent tie breaking.
    assert m['partial_ap_r0.1'] == pytest.approx(.5)
    assert m == ranking_metrics(y[::-1], p[::-1])
    assert ranking_metrics(np.ones(2), np.ones(2))['roc_auc'] is None


def test_pair_formula_is_bounded_and_monotone():
    pytest.importorskip('scipy')
    x = np.array([[0, 0], [0, 1], [1, 0], [1, 1.]])
    model = fit_pair(x, np.array([0, 0, 1, 1]), interaction=True)
    p = predict_pair(model, x)
    assert np.all((p >= 0) & (p <= 1))
    assert p[0] <= min(p[1:3]) <= max(p[1:3]) <= p[3]
    assert np.all(np.asarray(model['weights']) >= 0)


def test_geometry_support_does_not_change_labels():
    pytest.importorskip('sklearn')
    from vdjdb.motifs.tcremp import fit_clusters
    q = pl.DataFrame({'species': ['HomoSapiens']*6, 'gene': ['TRA']*6,
        'antigen.epitope': ['PEP']*6, 'junction_aa': ['CAAF','CACF','CADF','CAEF','CAGF','CAHF'],
        'v_call': ['TRAV1*01']*6, 'j_call': ['TRAJ1*01']*6,
        'duplicate_count': [1]*6, 'clonotype_id': [str(i) for i in range(6)]})
    x = np.array([[0,0],[.1,0],[.2,0],[.3,0],[.4,0],[3,0]])
    ref = fit_clusters(q, x, .5)
    before = ref.assign(q,x)
    got = ref.support(q,x)
    assert before.equals(ref.assign(q,x))
    assert got['core_support'][0] > got['core_support'][-1]
    assert got['core_density'][0] > got['core_density'][-1]


def test_optional_integration_preserves_default_and_hard_assignments():
    pytest.importorskip('pgmpy')
    from vdjdb.score.evidence import separated_motifs
    from vdjdb.score.pipeline import load_model, record_scores
    r = pl.DataFrame({'record_id':['r1','r2'], 'species':['HomoSapiens']*2,
        'antigen.epitope':['PEP']*2, 'mhc.a':['HLA-A*02:01']*2, 'mhc.b':['B2M']*2,
        'mhc.class':['MHCI']*2, 'reference.id':['P1','P2'], 'meta.subject.id':['',''],
        'vdjdb.score':[1,1]})
    c = pl.DataFrame({'record_id':['r1','r2'], 'gene':['TRA']*2, 'cdr3':['CAAF','CACF'],
        'v.segm':['TRAV1*01']*2, 'j.segm':['TRAJ1*01']*2})
    e = pl.DataFrame({'record_id':['r1','r2'], 'gene':['TRA']*2,
        'evidence_type':['motif_tcrnet','motif_tcrnet_nonmember'],
        'evidence_value':['',''], 'evidence_source':['',''], 'evidence_score':[1.,0.]})
    from vdjdb.score.evidence import integrate
    hard=integrate(r,c,e)
    near=separated_motifs(hard,r,c,e)
    assert near['a_tcrnet_hard'].to_list()==hard['a_tcrnet'].to_list()
    assert near['a_tcrnet'].to_list()==pytest.approx([5/6,1/6])
    assert near['b_tcrnet'].null_count()==2
    model,digest=load_model()
    default=record_scores(r,c,e,model,model_sha256=digest)
    explicit=record_scores(r,c,e,model,model_sha256=digest,continuous_membership=False)
    assert default.equals(explicit)
    option=record_scores(r,c,e,model,model_sha256=digest,continuous_membership=True)
    assert option['continuous_membership'].all()
    assert option.height==default.height
    with pytest.raises(ValueError,match='Choose'):
        record_scores(r,c,e,model,model_sha256=digest,continuous_membership=True,soft_membership=True)


def test_fractional_index_keeps_hard_pair_support():
    pytest.importorskip('pgmpy')
    from vdjdb.score.bayes import score_records
    from vdjdb.score.pipeline import load_model
    model,_=load_model()
    frame=pl.DataFrame({'a_tcrnet':[.8], 'b_tcrnet':[.7],
                        'a_tcrnet_hard':[1], 'b_tcrnet_hard':[1]})
    got=score_records(model,frame,soft_membership=True)[0]
    assert got['pair_supported']


def test_threshold_policy_falls_back_for_unseen_or_one_class_peptides():
    from vdjdb.validate.membership import apply_thresholds, fit_thresholds
    y=np.array([0,1,0,1,1])
    score=np.array([.1,.8,.2,.9,.7])
    peptide=np.array(['P','P','P','P','Q'])
    policy=fit_thresholds(y,score,peptide)
    assert 'P' in policy['epitope_thresholds'] and 'Q' not in policy['epitope_thresholds']
    got=apply_thresholds(policy,np.array(['P','Q','UNSEEN']))
    assert got[1]==got[2]==policy['global_threshold']


def test_geometry_support_preserves_missing_peptide_eligibility():
    pytest.importorskip('sklearn')
    from vdjdb.motifs.tcremp import fit_clusters
    q=pl.DataFrame({'species':['HomoSapiens']*5,'gene':['TRA']*5,
        'antigen.epitope':['PEP']*5,'junction_aa':['CAAF','CACF','CADF','CAEF','CAGF'],
        'v_call':['TRAV1*01']*5,'j_call':['TRAJ1*01']*5,
        'duplicate_count':[1]*5,'clonotype_id':[str(i) for i in range(5)]})
    x=np.c_[np.arange(5)/10,np.zeros(5)]
    ref=fit_clusters(q,x,.5)
    query=pl.concat([q.head(1),q.head(1).with_columns(pl.lit(None,dtype=pl.String).alias('antigen.epitope'))])
    got=ref.support(query,x[:2])
    assert got['eligible'].to_list()==[True,False]
    assert got['core_support'][-1] is None
