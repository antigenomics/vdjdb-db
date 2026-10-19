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


def _motifs(records: pl.DataFrame, chains: pl.DataFrame, evidence: pl.DataFrame) -> pl.DataFrame:
    return evidence.filter(pl.col('evidence_type').str.starts_with('motif_')).join(
        chains.select('record_id', 'gene', 'cdr3', 'v.segm', 'j.segm'),
        on=['record_id', 'gene'], validate='m:1').join(
        records.select('record_id', 'species', 'antigen.epitope'), on='record_id', validate='m:1')


def canonical(frame: pl.DataFrame) -> pl.DataFrame:
    """Gene-level receptor:pMHC keys; both junction anchors are retained."""
    if missing := sorted(set(KEY) - set(frame.columns)):
        raise ValueError(f'Interaction identity needs {missing}')
    frame = frame.with_columns(pl.col(KEY).cast(pl.String).fill_null('').str.strip_chars())
    if frame.select(pl.any_horizontal(pl.col(KEY).str.contains(SEP, literal=True)).any()).item():
        raise ValueError('Identity contains a reserved separator')
    return frame.with_columns(pl.col('a_v', 'a_j', 'b_v', 'b_j').str.split('*').list.first()).with_columns(
        pl.concat_str(KEY, separator=SEP).alias('interaction_key'))


def interactions(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """Attach canonical paired identities to records, retaining both junction anchors."""
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
    return canonical(full)


def integrate(records: pl.DataFrame, chains: pl.DataFrame, evidence: pl.DataFrame,
              queries: pl.DataFrame | None = None) -> pl.DataFrame:
    """Summarize completed build tables, without reinterpreting their annotations.

    Publication totals use the manuscript's chain/peptide/V-gene key. Donor totals count
    distinct nonempty publication/donor contexts for the integrated pair:pMHC. Motif
    flags mean membership exists among its gene-level build entries. Only explicit
    producer nonmembers yield zero; absent evidence remains unknown.
    """
    full = interactions(records, chains)
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
    motif = _motifs(records, chains, evidence)
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


def soft_motifs(features: pl.DataFrame, records: pl.DataFrame, chains: pl.DataFrame,
                evidence: pl.DataFrame) -> pl.DataFrame:
    """Optional Hamming-1 annotated-neighbour fractions, within each producer's scope.

    One batched native search uses substitutions only, including exact matches. Distinct
    annotated neighbours contribute once; unannotated neighbourhoods stay missing.
    Counts accompany each fraction. These are not calibrated binding probabilities.
    """
    import seqtree

    ref = _motifs(records, chains, evidence).with_columns(
        pl.col('v.segm', 'j.segm').str.split('*').list.first())
    queries = pl.concat([features.select('interaction_key', 'species', 'antigen.epitope',
        pl.lit(gene).alias('gene'), pl.col(prefix + '_junction').alias('query'),
        pl.col(prefix + '_v').alias('v.segm'), pl.col(prefix + '_j').alias('j.segm'))
        for prefix, gene in [('a', 'TRA'), ('b', 'TRB')]]).filter(pl.col('query') != '')
    refs = ref['cdr3'].unique().sort().to_list()
    qs = queries['query'].unique().sort().to_list()
    hits = (seqtree.Index.build(refs).search_batch(qs, seqtree.SearchParams(
        max_subs=1, max_ins=0, max_dels=0, max_total_edits=1), threads=1) if refs and qs else [])
    edges = pl.DataFrame([(q, refs[h.ref_id]) for q, group in zip(qs, hits, strict=bool(hits))
                          for h in group], schema={'query': pl.String, 'cdr3': pl.String}, orient='row')
    output = features
    for method in ['tcrnet', 'tcremp']:
        key = ['species', 'gene', 'antigen.epitope', 'cdr3']
        if method == 'tcremp':
            key += ['v.segm', 'j.segm']
        annotated = ref.filter(pl.col('evidence_type').is_in(
            ['motif_' + method, 'motif_' + method + '_nonmember'])).group_by(key).agg(
                (pl.col('evidence_type') == 'motif_' + method).any().alias('member'))
        totals = queries.join(edges, on='query', validate='m:m').join(
            annotated, on=key, how='inner', validate='m:1').group_by('interaction_key', 'gene').agg(
                pl.len().alias('neighbours'), pl.col('member').sum().alias('members'))
        for prefix, gene in [('a', 'TRA'), ('b', 'TRB')]:
            name = prefix + '_' + method
            part = totals.filter(pl.col('gene') == gene).select('interaction_key',
                (pl.col('members') / pl.col('neighbours')).alias(name),
                pl.col('members').alias(name + '_members'),
                pl.col('neighbours').alias(name + '_neighbours'))
            output = output.drop(name).join(part, on='interaction_key', how='left', validate='1:1')
    return output.sort('interaction_key')


def separated_motifs(features: pl.DataFrame, records: pl.DataFrame, chains: pl.DataFrame,
                     evidence: pl.DataFrame) -> pl.DataFrame:
    """Optional hard-order-preserving neighbour support for improved scoring.

    Keep original assignments in *_hard. Values are experimental interpolation
    weights for the existing soft inference, not calibrated membership probabilities.
    Only the original producer-eligible rows receive a continuous index.
    """
    names = [f'{p}_{m}' for p in ('a', 'b') for m in ('tcrnet', 'tcremp')]
    categorical(features.select(names))  # The preserved assignments must be binary.
    hard = features.select('interaction_key', *(pl.col(n).alias(n + '_hard') for n in names))
    near = soft_motifs(features, records, chains, evidence).join(
        hard, on='interaction_key', validate='1:1')
    return near.with_columns(((2*pl.col(n + '_hard') + pl.col(n))/3).alias(n) for n in names)


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
