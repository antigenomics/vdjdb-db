"""Optional pgmpy contract, run in the pinned confidence environment."""
import json

import pytest

pytest.importorskip('pgmpy')
from vdjdb.score.bayes import fit, predict, predict_marginals


def test_exact_missing_inference_roundtrip_and_deterministic_tree():
    rows = [{'a_studies': '2' if i % 2 else '0', 'b_tcrnet': '1' if i % 3 else '0'}
            for i in range(24)]
    model = fit(rows, [i % 2 for i in range(24)])
    restored = json.loads(json.dumps(model))
    p = predict(restored, [{}, {'a_studies': None}, {'a_studies': '2'}, {'a_studies': '0'}])
    assert p[0] == pytest.approx(.5)
    assert p[1] == p[0] and p[2] > p[0] > p[3]
    assert model == fit(rows, [i % 2 for i in range(24)])
    views = predict_marginals(model, [rows[1], {}])
    assert views[0]['TRA'] == pytest.approx(predict(model, [{'a_studies': rows[1]['a_studies']}])[0])
    assert views[0]['TRB'] == pytest.approx(predict(model, [{'b_tcrnet': rows[1]['b_tcrnet']}])[0])
    assert views[1] == {'paired': p[0], 'TRA': p[0], 'TRB': p[0]}
    with pytest.raises(ValueError):
        predict(model, [{'padj': '0'}])
    with pytest.raises(ValueError):
        predict(model, [{'iptm': 'ultra'}])


def test_unobserved_legal_evidence_states_do_not_create_support():
    rows = [{'a_studies': '2' if i % 2 else '0'} for i in range(24)]
    model = fit(rows, [int(i % 4 == 0) for i in range(24)])
    assert predict(model, [{'iptm': 'high', 'donors': '3+'}])[0] == pytest.approx(predict(model, [{}])[0])


def test_naive_network_keeps_exact_missing_inference_with_fewer_edges():
    rows = [{'a_studies': str(i % 2), 'b_tcrnet': str(int(i % 3 == 0))}
            for i in range(36)]
    model = fit(rows, [i % 2 for i in range(36)], structure='naive')
    assert len(model['edges']) == 9
    assert all(parent == 'validated' for parent, _ in model['edges'])
    assert predict(model, [{}, {'a_studies': None}]) == pytest.approx([.5, .5])
    assert predict(model, [{'a_studies': '1'}])[0] > .5
    with pytest.raises(ValueError):
        fit(rows, [i % 2 for i in range(36)], structure='unknown')


def test_unpredicted_structure_is_marginalized_in_observed_chain_context():
    rows = [{'a_studies': '1', 'iptm': 'high' if i % 2 else 'low'} for i in range(24)]
    for structure in ['tan', 'naive']:
        model = fit(rows, [i % 2 for i in range(24)], structure=structure)
        scores = predict(model, [
            {'a_studies': '1'}, {'a_studies': '1', 'iptm': None},
            {'a_studies': '1', 'iptm': '?'}, {'a_studies': '1', 'iptm': 'high'},
        ])
        assert scores[:3] == pytest.approx([scores[0]] * 3)
        assert scores[3] > scores[0]


def test_leaf_smoothing_preserves_predictions_without_leaf_evidence():
    from vdjdb.score.evidence import FEATURES

    rows = [{'a_studies': str(i % 2), 'iptm': 'high' if i % 3 else 'low'}
            for i in range(36)]
    labels = [i % 2 for i in range(36)]
    sizes = {name: 10. for name in [*FEATURES, 'validated']}
    sizes['iptm'] = 50.
    baseline = fit(rows, labels, structure='naive')
    refined = fit(rows, labels, structure='naive', equivalent_sample_size=sizes)
    queries = [{}, {'a_studies': '0'}, {'a_studies': '1'}]
    assert predict(refined, queries) == pytest.approx(predict(baseline, queries))
    assert refined['edges'] == baseline['edges']
    for invalid in [{'iptm': 50.}, {**sizes, 'iptm': float('nan')}]:
        with pytest.raises(ValueError):
            fit(rows, labels, equivalent_sample_size=invalid)


def test_structural_leaf_keeps_backbone_independent_of_grading():
    rows = [{'a_studies': str(i % 2), 'b_tcrnet': str(int(i % 3 == 0)),
             'iptm': 'high' if i % 4 else 'low'} for i in range(48)]
    labels = [i % 2 for i in range(48)]
    original = fit(rows, labels, structure='tan_leaf')
    regraded = fit([{**row, 'iptm': 'mid' if i % 3 else 'low'}
                    for i, row in enumerate(rows)], labels, structure='tan_leaf')
    assert len(original['edges']) == len(regraded['edges']) <= 17
    assert not any(parent == 'iptm' for parent, _ in original['edges'])
    assert [e for e in original['edges'] if e[1] != 'iptm'] == [
        e for e in regraded['edges'] if e[1] != 'iptm']
    assert predict(original, [{'a_studies': '1'}, {'b_tcrnet': '0'}]) == pytest.approx(
        predict(regraded, [{'a_studies': '1'}, {'b_tcrnet': '0'}]))


def test_independent_structure_leaf_has_only_class_parent():
    rows = [{'a_studies': str(i % 2), 'iptm': 'high' if i % 3 else 'low'}
            for i in range(36)]
    model = fit(rows, [i % 2 for i in range(36)], structure='tan_independent_iptm')
    assert [e for e in model['edges'] if 'iptm' in e] == [['validated', 'iptm']]
    assert model == fit(rows, [i % 2 for i in range(36)], structure='tan_independent_iptm')


def test_numeric_scoring_uses_saved_bounds_and_missing_columns():
    import polars as pl

    from vdjdb.score.bayes import predict_evidence, predict_marginals

    rows = [{'iptm': 'high' if i % 2 else 'mid'} for i in range(24)]
    model = fit(rows, [i % 2 for i in range(24)], structure='tan_independent_iptm')
    model['iptm_bounds'] = [.7, .9]
    actual = predict_evidence(model, pl.DataFrame({'iptm': [.85, None], 'id': ['x', 'y']}))
    assert actual == predict_marginals(model, [{'iptm': 'mid'}, {}])
    assert predict_evidence(model, pl.DataFrame({'id': ['z']})) == predict_marginals(model, [{}])
    assert predict_evidence(model, pl.DataFrame({'iptm': []})) == []


def test_soft_motif_membership_is_a_mixture_with_unchanged_hard_endpoints(monkeypatch):
    import polars as pl

    from vdjdb.score.bayes import predict_evidence

    rows = [{'a_tcrnet': str(i % 2), 'b_tcremp': str(int(i % 3 == 0))}
            for i in range(36)]
    model = fit(rows, [i % 2 for i in range(36)], structure='naive')
    hard = predict(model, [{'a_tcrnet': '0'}, {'a_tcrnet': '1'}])
    with pytest.raises(ValueError):
        predict(model, [{'a_tcrnet': .25}])
    assert predict(model, [{'a_tcrnet': .25}], soft_membership=True)[0] == pytest.approx(
        .75 * hard[0] + .25 * hard[1])
    assert predict(model, [{'a_tcrnet': 0.}, {'a_tcrnet': 1.}]) == pytest.approx(hard)
    both = predict(model, [{'a_tcrnet': str(a), 'b_tcremp': str(b)}
                           for a in (0, 1) for b in (0, 1)])
    assert predict(model, [{'a_tcrnet': .25, 'b_tcremp': .6}], soft_membership=True)[0] == pytest.approx(
        sum(w * p for w, p in zip([.3, .45, .1, .15], both, strict=True)))
    marginal = predict_marginals(model, [{'a_tcrnet': .25, 'b_tcremp': .6}], soft_membership=True)[0]
    assert marginal['TRA'] == pytest.approx(.75 * hard[0] + .25 * hard[1])
    with pytest.raises(ValueError):
        predict_evidence(model, pl.DataFrame({'a_tcrnet': [.25]}))
    numeric = predict_evidence(model, pl.DataFrame({'a_tcrnet': [.25]}), soft_membership=True)
    assert numeric[0]['paired'] == pytest.approx(.75 * hard[0] + .25 * hard[1])
    from vdjdb.score import bayes

    batch = [{'a_tcrnet': (i + 1) / 40, 'b_tcremp': .6} for i in range(30)]
    expected = predict_marginals(model, batch, soft_membership=True)
    monkeypatch.setattr(bayes, '_ASSIGNMENT_BATCH', 16)
    actual = predict_marginals(model, batch, soft_membership=True)
    for before, after in zip(expected, actual, strict=True):
        assert after == pytest.approx(before)
    for value in [-.1, 1.1, float('nan'), float('inf')]:
        with pytest.raises(ValueError):
            predict(model, [{'a_tcrnet': value}])


def test_record_score_retains_forecast_and_describes_support_without_changing_marginals():
    import polars as pl

    from vdjdb.score.bayes import predict_evidence, score_records

    rows = [{'a_tcrnet': str(i % 2)} for i in range(36)]
    model = fit(rows, [i % 2 for i in range(36)], structure='naive')
    evidence = pl.DataFrame({'a_tcrnet': [1] * 6, 'b_studies': [0, 1, 0, 0, None, 0],
                             'iptm': [None, None, .8, None, None, None],
                             'legacy_score': [None, None, None, 0, None, None],
                             'b_tcremp': [None, None, None, None, None, 1]})
    old = predict_evidence(model, evidence)
    scored = score_records(model, evidence)
    assert all(z['paired'] > .5 for z in old)
    assert [z['pair_supported'] for z in scored] == [False, True, True, True, False, True]
    for before, after in zip(old, scored, strict=True):
        assert {n: after[n] for n in before} == before
        assert after['record_score'] == before['paired']
        assert after['paired_lower'] <= after['record_score'] <= after['paired_upper']
        assert after['positive'] == (True if after['paired_lower'] >= .5 else
                                     False if after['paired_upper'] < .5 else None)
    unknown = score_records(model, pl.DataFrame({'id': ['x', 'y']}))
    assert len(unknown) == 2 and all(z['record_score'] == z['paired'] and
                                   z['positive'] is None for z in unknown)
    assert score_records(model, pl.DataFrame({'id': []})) == []


def test_missing_partner_retains_prediction_and_exposes_cross_threshold_range():
    import polars as pl

    from vdjdb.score.bayes import score_records

    rows, labels = [], []
    for rep in range(30):
        for a in (0, 1):
            for b in (0, 1):
                y = int(a + b + (rep % 3 == 0) >= 2)
                rows.append({'a_tcrnet': str(a), 'b_tcrnet': str(b),
                             'iptm': 'high' if y else 'low'})
                labels.append(y)
    model = fit(rows, labels, structure='naive')
    evidence = pl.DataFrame({'a_tcrnet': [1, 1], 'b_tcrnet': [None, 1],
                             'iptm': [None, .9]})
    z = score_records(model, evidence)
    assert z[0]['record_score'] == z[0]['paired'] > .5
    assert z[0]['paired_lower'] < .5 < z[0]['paired_upper']
    assert z[0]['positive'] is None
    assert z[1]['paired_lower'] == z[1]['paired'] == z[1]['paired_upper']
    assert z[1]['positive'] is True
    assert not z[0]['direct_pair_evidence'] and z[1]['direct_pair_evidence']
    assert z[0]['paired_lower'] <= min(z[0]['TRA'], z[0]['TRB'])
    assert z[0]['paired_upper'] >= max(z[0]['TRA'], z[0]['TRB'])


def test_missing_evidence_batches_reuse_queries_only_within_the_call(monkeypatch):
    from collections import Counter

    from pgmpy.inference import VariableElimination

    from vdjdb.score import bayes
    from vdjdb.score.pipeline import load_model

    artifact, _ = load_model(legacy=True)
    rows = [{'a_studies': None, 'b_studies': '1', 'donors': '1',
             'a_tcrnet': q, 'b_tcrnet': 1-q, 'a_tcremp': q/2,
             'b_tcremp': 1-q/2, 'iptm': None, 'legacy_score': '0'}
            for q in (.25, .5, .75)]
    expected = bayes._paired_ranges(artifact, rows)
    original_init = VariableElimination.__init__
    original_query = VariableElimination.query
    counts = Counter()
    patterns = Counter()

    def init(self, *args, **kwargs):
        counts['models'] += 1
        return original_init(self, *args, **kwargs)

    def query(self, variables, evidence=None, **kwargs):
        patterns[tuple(sorted((evidence or {}).items()))] += 1
        return original_query(self, variables, evidence=evidence, **kwargs)

    monkeypatch.setattr(VariableElimination, '__init__', init)
    monkeypatch.setattr(VariableElimination, 'query', query)
    monkeypatch.setattr(bayes, '_ASSIGNMENT_BATCH', 32)
    for _ in range(2):
        counts.clear()
        patterns.clear()
        assert bayes._paired_ranges(artifact, rows) == expected
        assert counts['models'] == 1
        assert patterns and max(patterns.values()) == 1


def test_frozen_inference_is_exact_across_hash_seeds():
    import os
    import subprocess
    import sys

    code = """import json
from vdjdb.score.bayes import predict
from vdjdb.score.pipeline import load_model
model, _ = load_model()
rows = [{'a_studies': '1', 'a_tcrnet': str(a), 'a_tcremp': str(b)}
        for a in (0, 1) for b in (0, 1)]
print(json.dumps(predict(model, rows)))
"""
    outputs = [subprocess.run([sys.executable, '-c', code], check=True, capture_output=True,
                             text=True, timeout=60,
                             env={**os.environ, 'PYTHONHASHSEED': seed}).stdout
               for seed in ('1', '4')]
    assert outputs[0] == outputs[1]
