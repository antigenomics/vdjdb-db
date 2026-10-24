"""Synthetic release checks, without private assay rows."""
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjdb.score.pipeline import load_model


def test_v32_default_and_both_legacy_versions():
    model, _ = load_model()
    assert model['score_version'] == '3.2'
    assert load_model(version='3.1')[0]['score_version'] == '3.1'
    assert load_model(legacy=True)[0].get('score_version', '3') == '3'
    with pytest.raises(ValueError):
        load_model(version='invalid')


def test_encoding_structure_and_bins():
    from vdjdb.score.v32 import encode, structure_policy
    records = pl.DataFrame({'record_id': ['r'], 'meta.structure.id': ['test']})
    links = pl.DataFrame({'record_id': ['r'], 'interaction_key': ['k']})
    evidence = pl.DataFrame(schema={'record_id': pl.String, 'evidence_type': pl.String})
    frame = pl.DataFrame({'interaction_key': ['k'], 'legacy_score': [3.], 'iptm': [.2]})
    transformed = structure_policy(frame, records, links, evidence)
    assert transformed['legacy_score'].to_list() == [2.]
    assert transformed['iptm'].to_list() == [1.]
    model = {'iptm_bounds': [.4, .8], 'motif_bounds': {'a_tcremp': [.2, .8]}}
    coded = encode(model, transformed.with_columns(pl.lit(.4).alias('a_tcremp')))
    assert coded['legacy_score'].to_list() == ['2+']
    assert coded['iptm'].to_list() == ['high']
    assert coded['a_tcremp'].to_list() == ['1']
    with pytest.raises(ValueError):
        encode(model, pl.DataFrame({'a_tcremp': [1.1]}))


def test_standalone_help_and_missing_private_path():
    script = Path('processing/fit_record_score_v3_2.py')
    env = {k: v for k, v in os.environ.items() if k != 'TCRVDB_PATH'}
    help_run = subprocess.run([sys.executable, str(script), '--help'],
                              env=env, capture_output=True, text=True)
    assert help_run.returncode == 0, help_run.stderr
    assert 'TCRVDB_PATH' in help_run.stdout
    missing = subprocess.run([sys.executable, str(script), '--tables', 'unused',
                              '--selected-records', 'unused',
                              '--selection-manifest', 'unused', '--out', 'unused'],
                             env=env, capture_output=True, text=True)
    assert missing.returncode != 0 and 'TCRVDB_PATH' in missing.stderr


def test_standalone_labels_and_independent_target():
    spec = importlib.util.spec_from_file_location('fit_v32', 'processing/fit_record_score_v3_2.py')
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    truth = pl.DataFrame({'tcr_id': ['a', 'b'], 'pos': [1, 0], 'ambiguous': [0, 1]})
    assert script.labels(truth)['label'].to_list() == [1, -1]
    with pytest.raises(ValueError):
        script.labels(truth.with_columns(pl.lit(2).alias('pos')))
    r = pl.DataFrame({'record_id': ['a', 'b', 'c'], 'pmhc_id': ['p'] * 3,
                      'reference.id': ['one', 'one', 'two']})
    c = pl.DataFrame({'record_id': ['a', 'b', 'c'], 'clonotype_id': ['x'] * 3, 'cdr3': ['CAAF'] * 3})
    assert script.targets(r, c)['label'].to_list() == [1, 1, 1]
    assert script.targets(r.head(2), c.head(2))['label'].to_list() == [0, 0]


def test_learned_scores_missing_context_and_interval():
    from vdjdb.score.v32 import score_records
    model, _ = load_model()
    frame = pl.DataFrame({'antigen.epitope': ['OTHER', 'GILGFVFTL'], 'legacy_score': [1, 2],
                          'a_studies': [None, 1], 'iptm': [.9, None]})
    rows = score_records(model, frame)
    assert len(rows) == 2
    assert all(0 <= r['paired_lower'] <= r['paired'] <= r['paired_upper'] <= 1 for r in rows)
    assert np.isfinite([r['TRA'] for r in rows]).all()
    with pytest.raises(ValueError, match='epitope'):
        score_records(model, frame.drop('antigen.epitope'))


def test_bootstrap_metric_oracle_and_alignment():
    from processing.audit_record_scores import audit, weighted_metrics
    from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_auc_score
    y, score = np.array([0, 1, 0, 1]), np.array([.2, .9, .9, .9])
    w = np.array([[1., 2., 1., 1.]])
    actual = weighted_metrics(y, score, w)[0]
    idx = np.repeat(np.arange(4), w[0].astype(int))
    precision, recall, _ = precision_recall_curve(y[idx], score[idx])
    assert actual == pytest.approx([average_precision_score(y[idx], score[idx]),
                                   auc(recall, precision), roc_auc_score(y[idx], score[idx])])
    a = pl.DataFrame({'endpoint': ['test'] * 4, 'epitope': ['E'] * 4, 'case': list('abcd'),
                      'label': y, 'score': score, 'predictor': ['old'] * 4})
    cases = pl.concat([a, a.with_columns(pl.lit('new').alias('predictor'))])
    _, differences = audit(cases, 'new', ['old'], draws=20)
    assert differences['difference'].to_list() == [0.] * 6
    with pytest.raises(ValueError, match='identical'):
        audit(cases.head(7), 'new', ['old'], draws=20)


def test_fixed_motif_join_uses_nonmissing_value_in_any_record_order(monkeypatch):
    from vdjdb.score import v32
    features = pl.DataFrame({'interaction_key': ['k'], 'legacy_score': [1.], 'iptm': [None],
                            **{n: [0.] for n in v32.MOTIFS}})
    links = pl.DataFrame({'record_id': ['r1', 'r2'], 'interaction_key': ['k', 'k']})
    monkeypatch.setattr(v32, 'interactions', lambda *args: links)
    monkeypatch.setattr(v32, 'integrate', lambda *args: features)
    records = links.select('record_id')
    evidence = pl.DataFrame(schema={'record_id': pl.String, 'evidence_type': pl.String})
    motif = pl.DataFrame({'record_id': ['r1', 'r2'], **{n: [None, .2] for n in v32.MOTIFS}})
    for source in [motif, motif.reverse()]:
        result = v32.fixed_features(records, None, evidence, source)
        assert result.select(v32.MOTIFS).row(0) == (.2,) * 4


@pytest.mark.parametrize('value', [4., 2.5, float('inf'), -1.])
def test_v32_rejects_malformed_confidence_before_collapse(value):
    from vdjdb.score.v32 import encode, structure_policy
    model, _ = load_model()
    frame = pl.DataFrame({'legacy_score': [value], 'iptm': [.5], 'interaction_key': ['k']})
    with pytest.raises(ValueError):
        encode(model, frame)
    records = pl.DataFrame({'record_id': ['r'], 'meta.structure.id': ['test']})
    links = records.select('record_id').with_columns(pl.lit('k').alias('interaction_key'))
    evidence = pl.DataFrame(schema={'record_id': pl.String, 'evidence_type': pl.String})
    with pytest.raises(ValueError):
        structure_policy(frame, records, links, evidence)


def test_no_predictor_evidence_keeps_every_epitope_context():
    from vdjdb.score.v32 import score_records
    model, _ = load_model()
    result = score_records(model, pl.DataFrame({'antigen.epitope': ['UNKNOWN', 'GILGFVFTL']}))
    assert len(result) == 2
    assert all(np.isfinite(r['record_score']) and r['positive'] is None for r in result)


def test_private_label_construct_mismatch_and_conflicting_flags():
    spec = importlib.util.spec_from_file_location('fit_v32', 'processing/fit_record_score_v3_2.py')
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    truth = pl.DataFrame({'tcr_id': ['x'], 'pos': [1], 'ambiguous': [0], 'epitope_aa': ['E']})
    selected = pl.DataFrame({'tcr_id': ['x'], 'epitope_aa': ['OTHER']})
    with pytest.raises(ValueError, match='construct'):
        script.labels(truth, selected)
    with pytest.raises(ValueError, match='exclusive'):
        script.labels(truth.with_columns(pl.lit(1).alias('ambiguous')))
