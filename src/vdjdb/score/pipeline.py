"""Produce record confidence v3 from completed tables and a frozen Bayesian model."""
from __future__ import annotations

import gzip
import hashlib
import json
from importlib.resources import files
from pathlib import Path

import polars as pl

from ..emit.vdjdb3 import read_table
from .bayes import score_records
from .evidence import KEY, integrate, interactions, producer_motifs, soft_motifs

SCORE_SCHEMA = {
    'paired': pl.Float64, 'TRA': pl.Float64, 'TRB': pl.Float64,
    'pair_supported': pl.Boolean, 'direct_pair_evidence': pl.Boolean,
    'record_score': pl.Float64, 'paired_lower': pl.Float64, 'paired_upper': pl.Float64,
    'positive': pl.Boolean,
}


def load_model(path: Path | None = None, *, legacy: bool = False,
               version: str | None = None) -> tuple[dict, str]:
    """Read frozen parameters; inference never fits or downloads a model."""
    if sum([path is not None, legacy, version is not None]) > 1:
        raise ValueError('Choose a custom model or legacy v3, not both')
    version = '3' if legacy else version or '3.2'
    names = {'3': 'record_v3.json', '3.1': 'record_v3_1.json',
             '3.2': 'record_v3_2.json', '3.2b': 'record_v3_2b.json'}
    if version not in names:
        raise ValueError('Unknown score version')
    name = names[version]
    raw = path.read_bytes() if path else files('vdjdb.score').joinpath(name).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def record_scores(records: pl.DataFrame, chains: pl.DataFrame, evidence: pl.DataFrame,
                  model: dict, *, model_sha256: str, soft_membership: bool = False,
                  continuous_membership: bool | None = None,
                  motif_scores: pl.DataFrame | None = None) -> pl.DataFrame:
    """Score distinct interactions once and join to every current record ID.

    Hard producer memberships are used. Identical gene-level receptor:pMHC records
    share integrated evidence, while every observational record keeps its own row.
    """
    if soft_membership and continuous_membership:
        raise ValueError('Choose soft_membership or continuous_membership, not both')
    if continuous_membership is None:
        continuous_membership = not soft_membership and model.get('membership_mode') == 'producer'
    learned = model.get('format') == 'vdjdb-learned-bn-1'
    if learned and model.get('score_version') != '3.2':
        raise ValueError('Use the separate v3.2b route for independent-study models')
    if learned and (soft_membership or continuous_membership is False):
        raise ValueError('v3.2a requires its fitted continuous motif evidence')
    if motif_scores is not None and not learned:
        raise ValueError('Fixed motif input is supported for v3.2a only')
    features = integrate(records, chains, evidence)
    if soft_membership:
        features = soft_motifs(features, records, chains, evidence)
    if continuous_membership and motif_scores is None:
        features = producer_motifs(features, records, chains, evidence,
                                   width=model.get('membership_width', 1/12))
    if learned:
        from .v32 import fixed_features, structure_policy
        from .v32 import score_records as learned_scores

        if motif_scores is not None:
            features = fixed_features(records, chains, evidence, motif_scores)
        else:
            features = structure_policy(features, records, interactions(records, chains), evidence)
        forecasts = learned_scores(model, features)
    else:
        forecasts = score_records(model, features,
                                  soft_membership=soft_membership or continuous_membership)
    scores = pl.DataFrame(forecasts, schema=SCORE_SCHEMA)
    values = features.drop(KEY).hstack(scores)
    ids = interactions(records, chains).select('record_id', 'interaction_key')
    result = ids.join(values, on='interaction_key', how='left', validate='m:1').drop(
        'interaction_key').with_columns(pl.lit(model.get('score_version', '3')).alias('score_version'),
        pl.lit(model_sha256).alias('model_sha256'),
        pl.lit(soft_membership or continuous_membership).alias('soft_membership')).sort('record_id')
    if learned:
        from .learning import predict_network
        from .v32 import encode, independent_features

        independent, independent_digest = load_model(version='3.2b')
        focal = independent_features(records, chains, evidence)
        focal = focal.with_columns(pl.Series('independent_validation_score',
            predict_network(independent, encode(independent, focal))))
        result = result.join(focal.select('record_id', 'independent_validation_score'),
                             on='record_id', validate='1:1').with_columns(
            pl.col('record_score').alias('validation_score'),
            pl.lit('3.2a').alias('validation_score_version'),
            pl.lit('3.2b').alias('independent_validation_score_version'),
            pl.lit(independent_digest).alias('independent_model_sha256'))
    if continuous_membership:
        result = result.with_columns(pl.lit(True).alias('continuous_membership'),
                                     pl.lit('separated-producer-v2').alias('motif_score_version'),
                                     pl.lit(model.get('membership_width', 1/12)).alias('motif_score_width'))
    if model.get('score_version') in ('3.1', '3.2') or model.get('decision_threshold', .5) != .5:
        result = result.with_columns(pl.lit(model.get('decision_threshold', .5)).alias('decision_threshold'))
    return result


def write_scores(frame: pl.DataFrame, out: Path) -> None:
    """Write probabilities at 12 significant figures, preserving full-precision decisions."""
    if frame['record_id'].null_count() or frame['record_id'].n_unique() != frame.height:
        raise ValueError('Score output requires unique nonmissing record IDs')
    out.mkdir(parents=True, exist_ok=True)
    probabilities = [name for name, kind in SCORE_SCHEMA.items() if kind == pl.Float64]
    probabilities += [name for name in ('validation_score', 'independent_validation_score')
                      if name in frame.columns]
    frame = frame.sort('record_id').with_columns(pl.col(probabilities).round_sig_figs(12))
    frame.write_parquet(out / 'record_scores.parquet')
    with (out / 'record_scores.tsv.gz').open('wb') as raw, gzip.GzipFile(
            filename='', fileobj=raw, mode='wb', mtime=0) as zipped:
        frame.write_csv(zipped, separator='\t', line_terminator='\n')


def run(tables: Path, out: Path, model_path: Path | None = None, *,
        soft_membership: bool = False, continuous_membership: bool | None = None,
        legacy: bool = False, threshold: float | None = None,
        version: str | None = None, motif_scores: Path | None = None) -> pl.DataFrame:
    """Recompute scores from this build's records, chains and completed motif evidence."""
    from threadpoolctl import threadpool_limits

    if motif_scores and motif_scores.resolve() == (out / 'record_scores.parquet').resolve():
        raise ValueError('Fixed motif input must not be overwritten by score output')
    model, digest = load_model(model_path, legacy=legacy, version=version)
    if threshold is not None:
        model = {**model, 'decision_threshold': threshold}
    records = read_table(tables, 'records')
    chains = read_table(tables, 'chains')
    evidence = read_table(tables, 'evidence')
    if not evidence['evidence_type'].str.starts_with('motif_').any():
        raise ValueError('Run vdjdb motifs on these tables before score-records')
    with threadpool_limits(limits=1):
        frame = record_scores(records, chains, evidence, model, model_sha256=digest,
                              soft_membership=soft_membership, continuous_membership=continuous_membership,
                              motif_scores=pl.read_parquet(motif_scores) if motif_scores else None)
    write_scores(frame, out)
    return frame
