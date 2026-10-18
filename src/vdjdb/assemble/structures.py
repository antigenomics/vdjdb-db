"""Optional structure links proved by complete identity or modeled variable domains."""
from __future__ import annotations

import os
import tarfile
from pathlib import Path

import polars as pl

from .master import HASH_FIELDS

CONTEXT = ('species', 'antigen.epitope', 'mhc.a', 'mhc.b', 'mhc.class')
IDENTITY = ('species', *HASH_FIELDS, 'mhc.class')


def _coordinates(path: Path, wanted: set[str]) -> pl.DataFrame:
    """Read the requested complexes in one archive pass, without extracting files."""
    parts = []
    seen = set()
    incomplete = 0
    with tarfile.open(path, mode='r|gz') as archive:
        for entry in archive:
            if not entry.isfile() or not entry.name.endswith('_aa_coordinates.tsv'):
                continue
            key = Path(entry.name).name.removesuffix('_aa_coordinates.tsv')
            if key not in wanted:
                continue
            if key in seen:
                raise ValueError('duplicate coordinate member for one structure hash')
            seen.add(key)
            stream = archive.extractfile(entry)
            assert stream is not None
            data = pl.read_csv(stream, separator='\t', columns=['Chain', 'Residue', 'ResNum'],
                               schema_overrides={'ResNum': pl.Int64})
            data = data.filter(pl.col('Chain').is_in(['TCR_alpha', 'TCR_beta', 'peptide']))
            if data.select('Chain', 'ResNum').n_unique() != data.height:
                raise ValueError('duplicate coordinate residue for one structure hash')
            if data.filter(~pl.col('Residue').fill_null('').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]$')
                           | pl.col('ResNum').is_null() | (pl.col('ResNum') < 1)).height:
                raise ValueError('coordinate residues require amino acids and positive residue numbers')
            bounds = data.group_by('Chain').agg(pl.len().alias('n'),
                pl.col('ResNum').min().alias('first'), pl.col('ResNum').max().alias('last'))
            if bounds.height != 3 or bounds.filter((pl.col('first') != 1) |
                    (pl.col('last') - pl.col('first') + 1 != pl.col('n'))).height:
                incomplete += 1
                continue
            parts.append(data.sort('Chain', 'ResNum').group_by('Chain', maintain_order=True)
                         .agg(pl.col('Residue').str.join('').alias('sequence'))
                         .with_columns(pl.lit(key).alias('tcr_pmhc_hash')))
    if incomplete:
        print(f"structure coordinates: skipped {incomplete:,} incomplete numbered members", flush=True)
    if not parts:
        return pl.DataFrame(schema={'tcr_pmhc_hash': pl.String, 'alpha_domain': pl.String,
                                    'beta_domain': pl.String, 'modeled_peptide': pl.String})
    out = pl.concat(parts).pivot(on='Chain', index='tcr_pmhc_hash', values='sequence').rename({
        'TCR_alpha': 'alpha_domain', 'TCR_beta': 'beta_domain', 'peptide': 'modeled_peptide'}, strict=False)
    return out.with_columns([pl.lit(None, pl.String).alias(c) for c in
                            ['alpha_domain', 'beta_domain', 'modeled_peptide'] if c not in out.columns])


def full_reference_root() -> Path:
    """Require a fixed local reference before any ARDA API or build can fetch one."""
    root = os.environ.get("ARDA_HOME")
    if not root:
        raise ValueError("full structure germline requires ARDA_HOME before process startup")
    database = Path(root).expanduser().resolve() / "database"
    if not (database / "vdj").is_dir():
        raise ValueError("full structure germline requires a local complete reference")
    from arda.paths import database_dir

    if database_dir().resolve() != database:
        raise ValueError("ARDA_HOME must be fixed before reference APIs are first used")
    return database


def _full_flanks(gene: str, organism: str) -> tuple[dict[str, str], dict[str, str]]:
    """Read the explicitly configured local ARDA reference through its published API."""
    base = full_reference_root() / "vdj" / organism
    if not all((base / name).is_file() for name in
               ("alleles.fasta", "markup.tsv", "cdr3_anchors.tsv")):
        raise ValueError(f"full structure germline requires a local complete {organism} reference")
    from vdjtools.model.reference import arda_full_germline

    full = arda_full_germline(gene, organism)
    return ({allele: seq[:anchor] for (segment, allele), (seq, anchor) in full.items() if segment == "V"},
            {allele: seq[anchor + 3:] for (segment, allele), (seq, anchor) in full.items() if segment == "J"})


def _domains(chains: pl.DataFrame, records: pl.DataFrame,
             full_germline: bool = False) -> pl.DataFrame:
    """Rebuild amino-acid domains using explicit alleles, never gene-level defaults."""
    from vdjtools.model import load_bundled, translate

    from ..annotate.contig import _flanks
    from ..annotate.junction import MODELS

    parts = []
    joined = chains.filter(pl.col('cdr3') != '').join(
        records.select('record_id', 'species'), on='record_id', validate='m:1')
    for (species, gene), part in sorted(joined.group_by('species', 'gene'), key=lambda item: item[0]):
        if species not in MODELS or gene not in ('TRA', 'TRB'):
            continue
        source, organism = MODELS[species]
        extra_pre, extra_suf = _full_flanks(gene, organism) if full_germline else ({}, {})
        pre, suf = _flanks(load_bundled(gene, source, organism=organism), gene)
        pre = extra_pre | pre
        suf = extra_suf | suf
        pre = {k: translate(v) for k, v in pre.items()}
        suf = {k: translate(v) for k, v in suf.items()}
        parts.append(part.select('record_id', 'gene', pl.concat_str(
            pl.col('v.segm').replace_strict(pre, default=None), pl.col('cdr3'),
            pl.col('j.segm').replace_strict(suf, default=None), ignore_nulls=False).alias('domain')))
    if not parts:
        return pl.DataFrame(schema={'record_id': pl.String, 'alpha_domain': pl.String,
                                    'beta_domain': pl.String})
    out = pl.concat(parts).pivot(on='gene', index='record_id', values='domain').rename({
        'TRA': 'alpha_domain', 'TRB': 'beta_domain'}, strict=False)
    return out.with_columns([pl.lit(None, pl.String).alias(c) for c in
                            ['alpha_domain', 'beta_domain'] if c not in out.columns])


def _protein_context(frame: pl.DataFrame) -> pl.DataFrame:
    # Later HLA fields refine synonymous coding or noncoding sequence; expression suffixes stay.
    return frame.with_columns([pl.col(c).str.replace(
        r'^(HLA-[A-Z0-9]+\*\d+:\d+)(?::\d+)*$', '${1}').alias(c) for c in ('mhc.a', 'mhc.b')])


def recover_links(records: pl.DataFrame, chains: pl.DataFrame, metadata: pl.DataFrame,
                  coordinates: Path, exact: pl.DataFrame, *, full_germline: bool = False) -> pl.DataFrame:
    """Additional source-hash/record links; require both chains and the recorded pMHC."""
    from ..curate.nomenclature import harmonise_segments

    if full_germline:
        full_reference_root()
    if not set(IDENTITY).issubset(metadata.columns):
        raise ValueError(f'coordinate recovery requires structure identity fields {list(IDENTITY)}')
    exact = exact.select('tcr_pmhc_hash', 'record_id')
    source = (metadata.select('tcr_pmhc_hash', *IDENTITY)
              .with_columns(pl.col(pl.String).fill_null('')).unique())
    if source.group_by('tcr_pmhc_hash').len().filter(pl.col('len') > 1).height:
        raise ValueError('one structure hash has incompatible metadata identities')
    source = source.join(exact.select('tcr_pmhc_hash').unique(), on='tcr_pmhc_hash', how='anti')
    if source.is_empty():
        return exact.sort('tcr_pmhc_hash', 'record_id')
    current = records.select('record_id', *CONTEXT)
    for gene, suffix in [('TRA', 'alpha'), ('TRB', 'beta')]:
        side = chains.filter(pl.col('gene') == gene).select('record_id',
            *[pl.col(a).alias(f'{b}.{suffix}') for a, b in
              [('cdr3', 'cdr3'), ('v.segm', 'v'), ('j.segm', 'j')]])
        current = current.join(side, on='record_id', how='inner', validate='1:1')
    renamed, _ = harmonise_segments(source)
    # A spelling synonym must not discard an allele suffix or resolve a multi-call to one gene.
    genes = ['v.alpha', 'j.alpha', 'v.beta', 'j.beta']
    originals = source.select('tcr_pmhc_hash', *[pl.col(c).alias(f'__source_{c}') for c in genes])
    renamed = renamed.join(originals, on='tcr_pmhc_hash', validate='1:1').filter(
        pl.all_horizontal([pl.col(c).str.extract(r'(\*.*)$', 1).fill_null('') ==
                           pl.col(f'__source_{c}').str.extract(r'(\*.*)$', 1).fill_null('')
                           for c in genes]))
    complete = renamed.filter(pl.all_horizontal([pl.col(c) != '' for c in IDENTITY]))
    aliases = complete.join(current, on=list(IDENTITY), how='inner').select('tcr_pmhc_hash', 'record_id')
    remaining = source.join(aliases.select('tcr_pmhc_hash').unique(), on='tcr_pmhc_hash', how='anti')
    if remaining.is_empty():
        return pl.concat([exact, aliases]).unique().sort('tcr_pmhc_hash', 'record_id')
    query = remaining.join(_coordinates(coordinates, set(remaining['tcr_pmhc_hash'])),
                           on='tcr_pmhc_hash', how='inner')
    query = _protein_context(query).filter((pl.col('modeled_peptide') == pl.col('antigen.epitope')) &
        pl.all_horizontal([pl.col(c) != '' for c in CONTEXT]))
    context = _protein_context(current.select('record_id', *CONTEXT))
    candidates = context.join(query.select(*CONTEXT).unique(), on=list(CONTEXT), how='semi')
    target = candidates.join(_domains(chains.join(candidates.select('record_id'), on='record_id', how='semi'),
                                     candidates, full_germline), on='record_id', how='inner')
    # Keep an exact beta domain as well as a modeled single extra terminal E, never strip all E's.
    variants = pl.concat([query, query.filter(pl.col('beta_domain').str.ends_with('E')).with_columns(
        pl.col('beta_domain').str.slice(0, pl.col('beta_domain').str.len_chars() - 1))])
    molecular = variants.drop_nulls(['alpha_domain', 'beta_domain']).join(target,
        on=[*CONTEXT, 'alpha_domain', 'beta_domain'], how='inner').select('tcr_pmhc_hash', 'record_id')
    links = pl.concat([exact, aliases, molecular]).unique().sort('tcr_pmhc_hash', 'record_id')
    print(f"structure recovery: {links['tcr_pmhc_hash'].n_unique() - exact['tcr_pmhc_hash'].n_unique():,} "
          'additional source hashes linked by complete identity or variable domains', flush=True)
    return links
