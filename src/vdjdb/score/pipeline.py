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
from .evidence import KEY, integrate, interactions, separated_motifs, soft_motifs

SCORE_SCHEMA = {
    'paired': pl.Float64, 'TRA': pl.Float64, 'TRB': pl.Float64,
    'pair_supported': pl.Boolean, 'direct_pair_evidence': pl.Boolean,
    'record_score': pl.Float64, 'paired_lower': pl.Float64, 'paired_upper': pl.Float64,
    'positive': pl.Boolean,
}


def load_model(path: Path | None = None) -> tuple[dict, str]:
    """Read frozen parameters; inference never fits or downloads a model."""
    raw = path.read_bytes() if path else files('vdjdb.score').joinpath('record_v3.json').read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def record_scores(records: pl.DataFrame, chains: pl.DataFrame, evidence: pl.DataFrame,
                  model: dict, *, model_sha256: str, soft_membership: bool = False,
                  continuous_membership: bool = False) -> pl.DataFrame:
    """Score distinct interactions once and join to every current record ID.

    Hard producer memberships are used. Identical gene-level receptor:pMHC records
    share integrated evidence, while every observational record keeps its own row.
    """
    if soft_membership and continuous_membership:
        raise ValueError('Choose soft_membership or continuous_membership, not both')
    features = integrate(records, chains, evidence)
    if soft_membership:
        features = soft_motifs(features, records, chains, evidence)
    if continuous_membership:
        features = separated_motifs(features, records, chains, evidence)
    scores = pl.DataFrame(score_records(model, features,
                          soft_membership=soft_membership or continuous_membership),
                          schema=SCORE_SCHEMA)
    values = features.drop(KEY).hstack(scores)
    ids = interactions(records, chains).select('record_id', 'interaction_key')
    result = ids.join(values, on='interaction_key', how='left', validate='m:1').drop(
        'interaction_key').with_columns(pl.lit('3').alias('score_version'),
        pl.lit(model_sha256).alias('model_sha256'),
        pl.lit(soft_membership or continuous_membership).alias('soft_membership')).sort('record_id')
    if continuous_membership:
        result = result.with_columns(pl.lit(True).alias('continuous_membership'),
                                     pl.lit('separated-neighbours-v1').alias('motif_score_version'))
    return result


def write_scores(frame: pl.DataFrame, out: Path) -> None:
    """Write probabilities at 12 significant figures, preserving full-precision decisions."""
    if frame['record_id'].null_count() or frame['record_id'].n_unique() != frame.height:
        raise ValueError('Score output requires unique nonmissing record IDs')
    out.mkdir(parents=True, exist_ok=True)
    probabilities = [name for name, kind in SCORE_SCHEMA.items() if kind == pl.Float64]
    frame = frame.sort('record_id').with_columns(pl.col(probabilities).round_sig_figs(12))
    frame.write_parquet(out / 'record_scores.parquet')
    with (out / 'record_scores.tsv.gz').open('wb') as raw, gzip.GzipFile(
            filename='', fileobj=raw, mode='wb', mtime=0) as zipped:
        frame.write_csv(zipped, separator='\t', line_terminator='\n')


def run(tables: Path, out: Path, model_path: Path | None = None, *,
        soft_membership: bool = False, continuous_membership: bool = False) -> pl.DataFrame:
    """Recompute scores from this build's records, chains and completed motif evidence."""
    from threadpoolctl import threadpool_limits

    model, digest = load_model(model_path)
    records = read_table(tables, 'records')
    chains = read_table(tables, 'chains')
    evidence = read_table(tables, 'evidence')
    if not evidence['evidence_type'].str.starts_with('motif_').any():
        raise ValueError('Run vdjdb motifs on these tables before score-records')
    with threadpool_limits(limits=1):
        frame = record_scores(records, chains, evidence, model, model_sha256=digest,
                              soft_membership=soft_membership, continuous_membership=continuous_membership)
    write_scores(frame, out)
    return frame
