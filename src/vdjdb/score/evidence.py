"""Compact confidence features from completed VDJdb build tables.

Use the shipped legacy score, publication/donor totals, motif membership and mean ipTM.
No chunk reconstruction, assay taxonomy, experiment tracing or allele-specific model.
"""
from __future__ import annotations

import math

import polars as pl

from ..identity.levels import SEP

CONTEXT = ['species', 'antigen.epitope', 'mhc.a', 'mhc.b', 'mhc.class']
RECEPTOR = [f'{c}_{f}' for c in ('a', 'b') for f in ('junction', 'v', 'j')]
KEY = CONTEXT + RECEPTOR
FEATURES = ['a_studies', 'b_studies', 'donors', 'a_tcrnet', 'b_tcrnet',
            'a_tcremp', 'b_tcremp', 'iptm', 'legacy_score']
COUNT_FEATURES = {'a_studies', 'b_studies', 'donors'}


def canonical(frame: pl.DataFrame) -> pl.DataFrame:
    """Gene-level receptor:pMHC keys; both junction anchors are retained."""
    if missing := sorted(set(KEY) - set(frame.columns)):
        raise ValueError(f'Interaction identity needs {missing}')
    frame = frame.with_columns(pl.col(KEY).cast(pl.String).fill_null('').str.strip_chars())
    if frame.select(pl.any_horizontal(pl.col(KEY).str.contains(SEP, literal=True)).any()).item():
        raise ValueError('Identity contains a reserved separator')
    return frame.with_columns(pl.col('a_v', 'a_j', 'b_v', 'b_j').str.split('*').list.first()).with_columns(
        pl.concat_str(KEY, separator=SEP).alias('interaction_key'))


def integrate(records: pl.DataFrame, chains: pl.DataFrame, evidence: pl.DataFrame,
              queries: pl.DataFrame | None = None) -> pl.DataFrame:
    """Summarize completed build tables, without reinterpreting their annotations.

    Publication totals use the manuscript's chain/peptide/V-gene key. Donor totals count
    distinct nonempty publication/donor contexts for the integrated pair:pMHC. Motif
    flags mean membership exists among its gene-level build entries. Only explicit
    producer nonmembers yield zero; absent evidence remains unknown.
    """
    if records['record_id'].n_unique() != records.height:
        raise ValueError('Repeated record_id')
    if chains.select('record_id', 'gene').unique().height != chains.height:
        raise ValueError('Repeated record/chain')
    full = records.select('record_id', *CONTEXT, 'reference.id', 'meta.subject.id', 'vdjdb.score')
    for prefix, gene in [('a', 'TRA'), ('b', 'TRB')]:
        part = chains.filter(pl.col('gene') == gene).select(
            'record_id', pl.col('cdr3').alias(prefix + '_junction'),
            pl.col('v.segm').alias(prefix + '_v'), pl.col('j.segm').alias(prefix + '_j'))
        full = full.join(part, on='record_id', how='left', validate='1:1')
    full = canonical(full)
    units = canonical(queries).select(*KEY, 'interaction_key') if queries is not None else full.select(
        *KEY, 'interaction_key').unique()
    if units['interaction_key'].n_unique() != units.height:
        raise ValueError('Query interactions must be unique')
    summary = full.group_by('interaction_key').agg(
        pl.col('vdjdb.score').cast(pl.Int64).max().alias('legacy_score'),
        pl.struct('reference.id', 'meta.subject.id').filter(
            pl.col('meta.subject.id').fill_null('') != '').n_unique().alias('donors'),
    ).with_columns(pl.when(pl.col('donors') > 0).then(pl.col('donors')).otherwise(None).alias('donors'))
    units = units.join(summary, on='interaction_key', how='left', validate='1:1')
    for prefix in ['a', 'b']:
        key = ['species', 'antigen.epitope', prefix + '_junction', prefix + '_v']
        support = full.filter(pl.col('reference.id').fill_null('') != '').group_by(key).agg(
            pl.col('reference.id').n_unique().alias(prefix + '_studies'))
        units = units.join(support, on=key, how='left', validate='m:1').with_columns(
            pl.when((pl.col(prefix + '_junction') != '') & (pl.col(prefix + '_v') != ''))
            .then(pl.col(prefix + '_studies').fill_null(0)).otherwise(None).alias(prefix + '_studies'))
    motif = evidence.filter(pl.col('evidence_type').str.starts_with('motif_')).join(
        chains.select('record_id', 'gene', 'cdr3', 'v.segm', 'j.segm'),
        on=['record_id', 'gene'], validate='m:1').join(
        records.select('record_id', 'species', 'antigen.epitope'), on='record_id', validate='m:1')
    for prefix, gene in [('a', 'TRA'), ('b', 'TRB')]:
        for method in ['tcrnet', 'tcremp']:
            key = ['species', 'antigen.epitope', prefix + '_junction']
            rename = {'cdr3': prefix + '_junction', 'v.segm': prefix + '_v', 'j.segm': prefix + '_j'}
            if method == 'tcremp':
                key += [prefix + '_v', prefix + '_j']
            part = motif.filter((pl.col('gene') == gene) & pl.col('evidence_type').is_in(
                ['motif_' + method, 'motif_' + method + '_nonmember'])).rename(rename).with_columns(
                    pl.col(prefix + '_v', prefix + '_j').str.split('*').list.first())
            name = prefix + '_' + method
            part = part.group_by(key).agg((pl.col('evidence_type') == 'motif_' + method)
                                          .any().cast(pl.Int64).alias(name))
            units = units.join(part, on=key, how='left', validate='m:1')
    models = evidence.filter((pl.col('evidence_type') == 'structure_model') &
                             (pl.col('evidence_value') == 'iptm')).join(
        full.select('record_id', 'interaction_key'), on='record_id', validate='m:1').select(
            'interaction_key', 'evidence_source', 'evidence_score').unique()
    if models.group_by('interaction_key', 'evidence_source').len().filter(pl.col('len') > 1).height:
        raise ValueError('Conflicting ipTM for one model')
    models = models.group_by('interaction_key').agg(pl.col('evidence_score').mean().alias('iptm'))
    return units.join(models, on='interaction_key', how='left', validate='1:1').sort('interaction_key')


def categorical(frame: pl.DataFrame, *, iptm_bounds: tuple[float, float] = (.5, .8)) -> pl.DataFrame:
    """Count bins 0/1/2/3+ and three ipTM grades at supplied bounds; '?' is unknown."""
    low, high = iptm_bounds
    if not all(math.isfinite(x) for x in iptm_bounds) or not 0 <= low < high <= 1:
        raise ValueError('ipTM bounds must be finite and ordered within zero to one')
    columns = []
    for name in frame.columns:
        value = pl.col(name).cast(pl.Float64, strict=True).fill_nan(None)
        if frame[name].cast(pl.Float64).is_infinite().any():
            raise ValueError(f'Infinite {name}')
        if name in COUNT_FEATURES:
            if frame.filter(value.is_not_null() & ((value < 0) | (value != value.floor()))).height:
                raise ValueError(f'Invalid count {name}')
            expr = pl.when(value >= 3).then(pl.lit('3+')).otherwise(value.cast(pl.Int64).cast(pl.String))
        elif name == 'iptm':
            if frame.filter(value.is_not_null() & ~value.is_between(0, 1)).height:
                raise ValueError('ipTM must be between zero and one')
            expr = pl.when(value < low).then(pl.lit('low')).when(value < high).then(
                pl.lit('mid')).otherwise(pl.lit('high'))
        else:
            allowed = [0, 1, 2, 3] if name == 'legacy_score' else [0, 1]
            if frame.filter(value.is_not_null() & ~value.is_in(allowed)).height:
                raise ValueError(f'Invalid state {name}')
            expr = value.cast(pl.Int64).cast(pl.String)
        columns.append(pl.when(value.is_null()).then(pl.lit('?')).otherwise(expr).alias(name))
    return frame.select(columns)
