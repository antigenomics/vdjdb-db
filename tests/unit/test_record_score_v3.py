"""Frozen v3 forecasts, record identity, compressed exports and matched comparisons."""
import gzip
from pathlib import Path

import polars as pl
import pytest

from vdjdb.score.pipeline import load_model, record_scores, write_scores
from vdjdb.validate.confidence import check_nonregression, metrics


def test_v31_default_and_explicit_legacy_model():
    model, digest = load_model()
    legacy, legacy_digest = load_model(legacy=True)
    assert model['score_version'] == '3.1'
    assert model['n_training'] == 1782 and len(model['training_epitopes']) == 5
    assert model['membership_mode'] == 'producer' and model['membership_width'] == 1/12
    assert model['decision_threshold'] == .26
    assert digest != legacy_digest and legacy['n_training'] == 614
    with pytest.raises(ValueError, match='Choose'):
        load_model(Path('unused.json'), legacy=True)


def test_decision_cutoff_changes_calls_without_changing_probabilities():
    pytest.importorskip('pgmpy')
    from vdjdb.score.bayes import score_records

    model, _ = load_model()
    frame = pl.DataFrame({'a_studies': [3], 'b_studies': [3], 'donors': [3],
        'a_tcrnet': [1], 'b_tcrnet': [1], 'a_tcremp': [1], 'b_tcremp': [1],
        'iptm': [.95], 'legacy_score': [3]})
    low = score_records({**model, 'decision_threshold': 0.}, frame)[0]
    high = score_records({**model, 'decision_threshold': 1.}, frame)[0]
    assert low['positive'] and not high['positive']
    assert {k: v for k, v in low.items() if k != 'positive'} == {
        k: v for k, v in high.items() if k != 'positive'}
    with pytest.raises(ValueError, match='threshold'):
        score_records({**model, 'decision_threshold': float('nan')}, frame)


def test_frozen_forecasts_and_missing_evidence():
    pytest.importorskip('pgmpy')
    from vdjdb.score.bayes import predict_evidence, score_records

    model, digest = load_model(legacy=True)
    assert digest == '10d70dad57f5e00153f315e138e135a3e33893435375bb2ddcc722148f129d53'
    assert model['n_training'] == 614 and model['structure'] == 'tan_independent_iptm'
    # Synthetic profiles, independent of assay identities or observed feature rows.
    frame = pl.DataFrame({'a_studies': [0, 1, 2, None, 3], 'b_studies': [0, 1, 1, None, 3],
                         'a_tcrnet': [0, 1, 1, None, 1], 'b_tcremp': [0, 0, 1, None, 1],
                         'iptm': [None, .4, .8, None, .95]})
    expected = [.2010713642247576, .16594253154766386, .9609077335708823,
                .4983974358974359, .9910220432098391]
    actual = score_records(model, frame)
    assert [r['record_score'] for r in actual] == pytest.approx(expected, abs=1e-14)
    assert [r['positive'] for r in actual] == [None, None, True, None, True]
    assert [{k: r[k] for k in ('paired', 'TRA', 'TRB')} for r in actual] == predict_evidence(model, frame)
    for row in actual:
        assert row['paired_lower'] <= row['record_score'] <= row['paired_upper']
    with pytest.raises(ValueError):
        score_records(model, pl.DataFrame({'a_tcrnet': [.25]}))


def test_record_join_and_deterministic_compressed_export(tmp_path):
    pytest.importorskip('pgmpy')
    records = pl.DataFrame({'record_id': ['r2', 'r1'], 'species': ['HomoSapiens'] * 2,
        'antigen.epitope': ['PEPTIDE'] * 2, 'mhc.a': ['HLA-A*02:01'] * 2, 'mhc.b': ['B2M'] * 2,
        'mhc.class': ['MHCI'] * 2, 'reference.id': ['P1', 'P2'], 'meta.subject.id': ['', ''],
        'vdjdb.score': [1, 1]})
    chains = pl.DataFrame({'record_id': ['r1', 'r2'], 'gene': ['TRB'] * 2, 'cdr3': ['CASSF'] * 2,
                          'v.segm': ['TRBV1*01'] * 2, 'j.segm': ['TRBJ1*01'] * 2})
    evidence = pl.DataFrame(schema={'record_id': pl.String, 'gene': pl.String,
        'evidence_type': pl.String, 'evidence_value': pl.String, 'evidence_score': pl.Float64,
        'evidence_source': pl.String})
    model, digest = load_model()
    scored = record_scores(records, chains, evidence, model, model_sha256=digest)
    assert scored['score_version'].unique().to_list() == ['3.1']
    assert scored['continuous_membership'].all()
    assert scored['motif_score_version'].unique().to_list() == ['separated-producer-v2']
    assert scored['record_id'].to_list() == ['r1', 'r2']
    assert scored['record_score'][0] == scored['record_score'][1]
    assert scored['b_studies'].to_list() == [2, 2]
    assert scored['a_studies'].null_count() == 2
    assert not scored['direct_pair_evidence'].any() and not scored['pair_supported'].any()
    write_scores(scored, tmp_path)
    before = (tmp_path / 'record_scores.tsv.gz').read_bytes()
    write_scores(scored.reverse(), tmp_path)
    assert (tmp_path / 'record_scores.tsv.gz').read_bytes() == before
    exported = pl.read_parquet(tmp_path / 'record_scores.parquet')
    probabilities = ['paired', 'TRA', 'TRB', 'record_score', 'paired_lower', 'paired_upper']
    assert exported.drop(probabilities).equals(scored.drop(probabilities))
    for name in probabilities:
        assert exported[name].to_list() == pytest.approx(scored[name].to_list(), rel=1e-11, abs=0)
    assert gzip.decompress(before).decode().splitlines()[0].startswith('record_id\t')
    with pytest.raises(ValueError):
        write_scores(pl.concat([scored, scored]), tmp_path)


def test_export_precision_ignores_machine_roundoff_and_preserves_small_probabilities(tmp_path):
    frame = pl.DataFrame({'record_id': ['r1'], 'paired': [.123456789012341],
        'record_score': [.123456789012341], 'TRA': [1.23456789012345e-18], 'TRB': [.4],
        'paired_lower': [.1], 'paired_upper': [.9], 'iptm': [.81234567890123456]})
    write_scores(frame, tmp_path)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    write_scores(frame.with_columns(pl.col('paired', 'record_score') + 2e-16), tmp_path)
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before
    exported = pl.read_parquet(tmp_path / 'record_scores.parquet')
    assert exported['paired'][0] == .123456789012
    assert exported['TRA'][0] == pytest.approx(1.23456789012e-18, rel=1e-15, abs=0)
    assert exported['iptm'][0] == frame['iptm'][0]


def test_binary_metrics_count_ties_and_reject_unmatched_cases():
    pytest.importorskip('sklearn')
    a = pl.DataFrame({'case': [1, 2, 3, 4], 'epitope': ['PEPTIDE'] * 4, 'fold': [0, 0, 1, 1],
        'y': [0, 0, 1, 1], 'score': [0., 1., 1., 1.], 'predicted': [0, 1, 1, 1],
        'predictor': ['motif'] * 4})
    table = metrics(pl.concat([a, a.with_columns(pl.lit('v3').alias('predictor'))]))
    check_nonregression(table)
    row = table.row(0, named=True)
    assert (row['tp'], row['fp'], row['fn'], row['tn']) == (2, 1, 0, 1)
    assert row['balanced_accuracy'] == .75 and row['fpr'] == .5
    assert row['ap'] == pytest.approx(2 / 3)
    assert row['false_discovery_fraction'] == pytest.approx(1 / 3)
    with pytest.raises(ValueError, match='identical cases'):
        metrics(pl.concat([a, a.with_columns(pl.lit('v3').alias('predictor')).head(3)]))
    with pytest.raises(ValueError, match='regresses'):
        check_nonregression(table.with_columns(pl.when(pl.col('predictor') == 'v3').then(.1)
                                              .otherwise(pl.col('ap')).alias('ap')))


def test_reviewed_aggregate_nonregression():
    table = pl.read_csv(Path(__file__).resolve().parents[2] /
                        'docs/_static/record-score-v3/metrics.tsv', separator='\t')
    check_nonregression(table)
    assert table.filter(pl.col('predictor') == 'v3')['n'].to_list() == [193, 421]
