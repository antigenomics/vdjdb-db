"""Legacy web co-deliverable, projected from shipped evidence without changing legacy exports."""
from pathlib import Path

import polars as pl

from .legacy import _CSV, build_default
from .vdjdb3 import joined
from ..io.text import write_csv, write_text
from ..schema import EVIDENCE_COLUMNS, VDJDB_WEB_COLUMNS
from ..schema.fields import render_meta


def build(tables: dict[str, pl.DataFrame]) -> pl.DataFrame:
    base = build_default(tables['records'], tables['chains'], identity=True)
    flags = joined(tables).with_columns(
        (pl.col('method.verification').fill_null('').str.strip_chars() != '')
        .alias('evidence.validation.same.study')).select('record_id', 'gene',
        'evidence.validation.same.study', 'evidence.validation.independent',
        'evidence.structure.native')
    # Model metrics are long evidence rows; native entries never supply predicted quality.
    ev = tables['evidence'].filter(pl.col('evidence_type') == 'structure_model')
    # Recovered model associations may use an older source hash. Keep an exact match;
    # otherwise require a unique model rather than pointing a badge at an arbitrary map.
    links = ev.select('record_id', pl.col('evidence_source').alias('__model')).unique()
    native = tables['evidence'].filter(pl.col('evidence_type') == 'structure_native').with_columns(
        pl.col('evidence_id').str.extract(r':([0-9a-f]{64})$', 1).alias('__model'))
    native = native.drop_nulls('__model').join(
        tables['records'].select('record_id', 'meta.structure.id'), on='record_id', validate='m:1')
    native = native.with_columns((pl.col('evidence_value').str.to_lowercase() ==
        pl.col('meta.structure.id').str.to_lowercase()).alias('__cited'))
    native = native.filter(pl.col('__cited') | ~pl.col('__cited').any().over('record_id'))
    native = native.select('record_id', '__model').unique()
    links = pl.concat([links.join(native.select('record_id').unique(), on='record_id', how='anti'), native])
    candidates = base.select('record_id', 'TCR_hash').unique().join(links, on='record_id')
    candidates = candidates.with_columns((pl.col('TCR_hash') == pl.col('__model')).alias('__exact'))
    candidates = candidates.filter(pl.col('__exact') | ~pl.col('__exact').any().over('record_id'))
    ambiguous = candidates.group_by('record_id').len().filter(pl.col('len') > 1)
    if ambiguous.height:
        from ..config import repo_root
        choices = pl.read_csv(repo_root() / 'rules/legacy_structure_representatives.tsv',
                              separator='\t').select('record_id', 'historical_hash', 'selected_hash')
        review = candidates.join(ambiguous.select('record_id'), on='record_id').join(
            choices, on='record_id', how='left', validate='m:1')
        selected = review.filter((pl.col('historical_hash') == pl.col('TCR_hash')) &
                                 (pl.col('selected_hash') == pl.col('__model')))
        if selected.height != ambiguous.height:
            raise ValueError('ambiguous structure associations require a current reviewed selection')
        candidates = pl.concat([
            candidates.join(ambiguous.select('record_id'), on='record_id', how='anti'),
            selected.select(candidates.columns)])
    base = base.join(candidates.select('record_id', '__model'), on='record_id', how='left',
                     validate='m:1').with_columns(
        pl.coalesce('__model', 'TCR_hash').alias('TCR_hash')).drop('__model')
    ev = ev.join(candidates.select('record_id', '__model'), on='record_id').filter(
        pl.col('evidence_source') == pl.col('__model'))
    metrics = ev.group_by('record_id').agg(
        ((pl.col('evidence_value') == 'num_contacts') & (pl.col('evidence_score') > 0))
        .any().alias('evidence.structure.contacts'),
        ((pl.col('evidence_value') == 'ranking_confidence') & (pl.col('evidence_score') >= .85))
        .any().alias('evidence.structure.quality'))
    return (base.join(flags, on=['record_id', 'gene'], how='left', validate='1:1')
        .join(metrics, on='record_id', how='left', validate='m:1')
        .with_columns(pl.col(*EVIDENCE_COLUMNS).fill_null(False))
        .select(VDJDB_WEB_COLUMNS))


def write_all(tables: dict[str, pl.DataFrame], out: Path, *, motifs: Path | None = None) -> dict[str, Path]:
    out.mkdir(parents=True, exist_ok=True)
    data = write_csv(build(tables), out / 'vdjdb.txt', **_CSV)
    meta = write_text(out / 'vdjdb.meta.txt', render_meta('vdjdb-web'))
    from ..config import repo_root
    import shutil
    shutil.copyfile(repo_root() / 'rules/legacy_native_representatives.tsv',
                    out / 'native_representatives.tsv')
    if motifs is not None:
        from .web_charts import write
        write(motifs, out / 'motif-charts')
    return {p.name: p for p in (data, meta)}
