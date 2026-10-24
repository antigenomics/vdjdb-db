"""v3.2a assay forecasts and v3.2b independent-study forecasts.

Fitting is offline in processing/. Continuous motif evidence comes from the
unchanged producer pipeline or an explicitly supplied fixed reference input.
"""
import math

import numpy as np
import polars as pl

from . import confidence
from .evidence import FEATURES, categorical, integrate, interactions
from .learning import predict_network

MOTIFS = [n for n in FEATURES if n.endswith(('tcrnet', 'tcremp'))]


def collapse_confidence(frame):
    """Validate the original 0-3 vocabulary before collapsing 2 and 3."""
    if 'legacy_score' not in frame.columns:
        return frame
    categorical(frame.select('legacy_score'))
    return frame.with_columns(pl.col('legacy_score').cast(pl.Float64).fill_nan(None).clip(upper_bound=2))


def structure_policy(features, records, links, evidence):
    """Collapse confidence and override interactions with experimental structures."""
    features = collapse_confidence(features)
    native = evidence.filter(pl.col('evidence_type') == 'structure_native').select('record_id')
    if 'meta.structure.id' in records.columns:
        native = pl.concat([native, records.filter(
            pl.col('meta.structure.id').fill_null('') != '').select('record_id')])
    keys = links.join(native.unique(), on='record_id', how='semi').select('interaction_key').unique()
    return features.join(keys.with_columns(pl.lit(True).alias('_native')),
                         on='interaction_key', how='left', validate='1:1').with_columns(
        pl.when(pl.col('_native').fill_null(False)).then(2.)
        .otherwise(pl.col('legacy_score').clip(upper_bound=2)).alias('legacy_score'),
        pl.when(pl.col('_native').fill_null(False)).then(1.)
        .otherwise(pl.col('iptm')).alias('iptm')).drop('_native')


def fixed_features(records, chains, evidence, motif_scores, units=None):
    """Join declared fixed motif values; keep unmatched interactions' hard evidence."""
    full = interactions(records, chains)
    features = integrate(records, chains, evidence, units)
    saved = motif_scores.select('record_id', *MOTIFS).join(
        full.select('record_id', 'interaction_key'), on='record_id', validate='1:1')
    conflicts = saved.group_by('interaction_key').agg(
        *(pl.col(n).drop_nulls().n_unique().alias(n) for n in MOTIFS))
    if conflicts.filter(pl.any_horizontal(pl.col(MOTIFS) > 1)).height:
        raise ValueError('Conflicting fixed motif values for one interaction')
    saved = saved.group_by('interaction_key').agg(
        *(pl.col(n).drop_nulls().first().alias(n + '_fixed') for n in MOTIFS))
    features = features.join(saved, on='interaction_key', how='left', validate='1:1').with_columns(
        *(pl.coalesce(n + '_fixed', n).alias(n) for n in MOTIFS)).drop([n + '_fixed' for n in MOTIFS])
    return structure_policy(features, records, full, evidence)


def encode(model, frame):
    """Apply the artifact's fixed bins; missing remains unobserved at inference."""
    motifs = model.get('motif_bounds', {})
    cols = [n for n in FEATURES if n in frame.columns and n not in motifs]
    frame = collapse_confidence(frame)
    result = categorical(frame.select(cols), iptm_bounds=model['iptm_bounds'])
    for name, (low, high) in motifs.items():
        if name not in frame.columns:
            continue
        v = pl.col(name).cast(pl.Float64).fill_nan(None)
        if frame.filter(v.is_not_null() & ~v.is_between(0, 1)).height:
            raise ValueError(f'Invalid motif value {name}')
        result = result.hstack(frame.select(
            pl.when(v.is_null()).then(pl.lit('?')).when(v <= low).then(pl.lit('0'))
            .when(v < .5).then(pl.lit('1')).when(v <= high).then(pl.lit('2'))
            .otherwise(pl.lit('3')).alias(name)))
    if 'legacy_score' in result.columns:
        result = result.with_columns(pl.col('legacy_score').replace({'2': '2+'}))
    if result.width == 0:
        result = pl.DataFrame({n: ['?'] * frame.height for n in model['states']})
    return result.select([n for n in FEATURES if n in result.columns])


def score_records(model, frame):
    """Paired and chain views, completion envelope, descriptive pair flags."""
    if 'antigen.epitope' not in frame.columns:
        raise ValueError('Aligned antigen.epitope context required')
    threshold = model.get('decision_threshold', .5)
    if not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('Decision threshold must be finite and in [0,1]')
    x = encode(model, frame)
    context = frame['antigen.epitope'].to_list()
    paired = np.asarray(predict_network(model, x, epitopes=context, bounds=True))
    a, b = [np.asarray(predict_network(model, x.with_columns(
        *(pl.lit('?').alias(n) for n in x.columns if not n.startswith(prefix))),
        epitopes=context)) for prefix in ('a_', 'b_')]

    def value(name):
        return pl.col(name) if name in frame.columns else pl.lit(None)

    has_pair = (pl.all_horizontal(pl.col('a_junction', 'b_junction').fill_null('') != '')
                if {'a_junction', 'b_junction'} <= set(frame.columns) else pl.lit(True))
    direct = has_pair & (value('legacy_score').is_not_null() | (value('donors') > 0) |
                         value('iptm').is_not_null())
    supported = [pl.any_horizontal(
        value(n + '_hard' if n + '_hard' in frame.columns else n) > 0 if n.endswith('studies')
        else value(n + '_hard' if n + '_hard' in frame.columns else n) >= .5
        for n in FEATURES if n.startswith(prefix)) for prefix in ('a_', 'b_')]
    flags = frame.with_columns(direct.fill_null(False).alias('direct'),
                         (has_pair & (direct | (supported[0] & supported[1])))
                         .fill_null(False).alias('support')).select('direct', 'support')
    rows = []
    for (p, low, high), aa, bb, direct, support in zip(
            paired, a, b, flags['direct'], flags['support'], strict=True):
        if not direct:
            low, high = min(low, aa, bb), max(high, aa, bb)
        rows.append(dict(paired=float(p), TRA=float(aa), TRB=float(bb), record_score=float(p),
                         paired_lower=float(low), paired_upper=float(high),
                         direct_pair_evidence=bool(direct), pair_supported=bool(support),
                         positive=True if low >= threshold else False if high < threshold else None))
    return rows


def independent_features(records, chains, evidence):
    """Focal-record evidence only; no counts or maxima across publications."""
    reads = pl.col('method.frequency.count').cast(pl.Int64, strict=False)
    f = records.with_columns(confidence.frequency().alias('__freq'),
        pl.coalesce(reads, confidence.cell_count()).alias('__count'), reads.alias('__reads'),
        pl.col('method.identification').str.strip_chars().str.to_lowercase().alias('__ident'),
        pl.col('method.verification').str.strip_chars().str.to_lowercase().alias('__verif'))
    f = f.select('record_id', 'antigen.epitope',
                 confidence.row_score().clip(upper_bound=2).alias('legacy_score'))
    native = pl.concat([records.filter(pl.col('meta.structure.id').fill_null('') != '').select('record_id'),
        evidence.filter(pl.col('evidence_type') == 'structure_native').select('record_id')]).unique()
    iptm = evidence.filter((pl.col('evidence_type') == 'structure_model') &
                          (pl.col('evidence_value') == 'iptm')).select(
                              'record_id', 'evidence_source', 'evidence_score').unique()
    if iptm.group_by('record_id', 'evidence_source').len().filter(pl.col('len') > 1).height:
        raise ValueError('Conflicting focal structure models')
    f = f.join(iptm.group_by('record_id').agg(pl.col('evidence_score').mean().alias('iptm')),
               on='record_id', how='left', validate='1:1').join(
                   native.with_columns(pl.lit(True).alias('_native')),
                   on='record_id', how='left', validate='1:1')
    f = f.with_columns(
        pl.when(pl.col('_native').fill_null(False)).then(1.).otherwise(pl.col('iptm')).alias('iptm'),
        pl.when(pl.col('_native').fill_null(False)).then(2).otherwise(pl.col('legacy_score')).alias('legacy_score'))
    for prefix, gene in [('a', 'TRA'), ('b', 'TRB')]:
        for method in ['tcrnet', 'tcremp']:
            part = evidence.filter((pl.col('gene') == gene) & pl.col('evidence_type').is_in(
                ['motif_' + method, 'motif_' + method + '_nonmember'])).group_by('record_id').agg(
                    (pl.col('evidence_type') == 'motif_' + method).any()
                    .cast(pl.Int64).alias(prefix + '_' + method))
            f = f.join(part, on='record_id', how='left', validate='1:1')
    return f.drop('_native')
