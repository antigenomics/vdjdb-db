"""Proofread current structure associations and deposited coordinate-file coverage.

The reviewed dataset mapping is an input, independent of automatic evidence assembly.
Archives are streamed without extraction. No molecular inference or source editing occurs.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import tarfile
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.request import urlopen

import polars as pl

from ..emit.vdjdb3 import read_table

REVISION = 'eb3dd781b8dbeaa4ce347e3761f7bdb1be44cb57'
BASE_URL = f'https://huggingface.co/datasets/isalgo/vdjdb_structure_models/resolve/{REVISION}/'


def _read(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator='\t', infer_schema_length=0).fill_null('')


def verify_inputs(root: Path, manifest: Path, *, fetch: bool = False,
                  names: tuple[str, ...] | None = None) -> int:
    """Verify pinned raw inputs; optionally stage missing files atomically from Hugging Face."""
    inputs = _read(manifest)
    if names is not None:
        if set(names) - set(inputs['path']):
            raise ValueError('Requested input is absent from the pinned manifest')
        inputs = inputs.filter(pl.col('path').is_in(names))
    for item in inputs.iter_rows(named=True):
        relative = Path(item['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Input path must be relative to the dataset root')
        path = root / relative
        if not path.is_file():
            if not fetch:
                raise ValueError(f'Missing structure input {relative}; use --fetch to stage pinned inputs')
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = None
            try:
                with NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                    temporary = Path(stream.name)
                    digest = hashlib.sha256()
                    with urlopen(BASE_URL + relative.as_posix(), timeout=60) as response:
                        while block := response.read(1024 * 1024):
                            digest.update(block)
                            stream.write(block)
                if digest.hexdigest() != item['sha256']:
                    raise ValueError(f'Download checksum mismatch: {relative}')
                temporary.replace(path)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        with path.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != item['sha256']:
                raise ValueError(f'Structure input checksum mismatch: {relative}')
    return inputs.height


def check_links(records: pl.DataFrame, chains: pl.DataFrame, links: pl.DataFrame,
                source: pl.DataFrame) -> None:
    """Fail on stale IDs, changed sequence/context snapshots or incorrect match flags."""
    if links.select('source_row', 'record_id').n_unique() != links.height:
        raise ValueError('Duplicate structure source/record link')
    current = records.select('record_id', *[
        pl.col(name).cast(pl.String).alias('current_' + name.replace('.', '_'))
        for name in ['species', 'antigen.epitope', 'mhc.a', 'mhc.b', 'mhc.class',
                     'reference.id', 'meta.structure.id', 'chunk.file', 'chunk.row']])
    for gene, side in [('TRA', 'alpha'), ('TRB', 'beta')]:
        current = current.join(chains.filter(pl.col('gene') == gene).select('record_id',
            pl.col('cdr3').alias('current_repaired_' + side),
            pl.col('cdr3.original').alias('current_original_' + side)),
            on='record_id', how='left', validate='1:1')
    current = current.with_columns(pl.col(pl.String).fill_null(''))
    checked = links.select('source_row', 'record_id').join(current, on='record_id',
        how='left', validate='m:1')
    if checked['current_species'].null_count():
        raise ValueError('Structure links contain retired or unknown record IDs')
    expected = links.select(checked.columns)
    if not checked.sort('source_row', 'record_id').equals(expected.sort('source_row', 'record_id')):
        raise ValueError('Structure mapping differs from current sequence or context; reconcile the dataset')
    source_fields = ['tcr_pmhc_hash', 'is_native', 'meta.structure.id']
    joined = links.join(source.select('source_row', 'cdr3.alpha', 'cdr3.beta', 'antigen.epitope',
        'species', 'mhc.a', 'mhc.b', 'mhc.class', *[
            pl.col(name).alias('source_' + name) for name in source_fields]),
        on='source_row', how='left', validate='m:1')
    if joined['cdr3.alpha'].null_count():
        raise ValueError('Unknown structure metadata source row')
    if joined.filter(pl.any_horizontal(
            pl.col(name) != pl.col('source_' + name) for name in source_fields)).height:
        raise ValueError('Structure link source identity differs from original metadata')
    for state in ['original', 'repaired']:
        match = ((pl.col('cdr3.alpha') == pl.col('current_' + state + '_alpha')) &
                 (pl.col('cdr3.beta') == pl.col('current_' + state + '_beta')) &
                 (pl.col('antigen.epitope') == pl.col('current_antigen_epitope')))
        if joined.filter(~pl.col('matches_' + state).is_in(['true', 'false']) |
                         (match != (pl.col('matches_' + state) == 'true'))).height:
            raise ValueError(f'Incorrect {state} structure match flag')
    def protein(name):
        return pl.col(name).str.replace(r'^H-2-?', 'H2-').str.replace(
            r'^(HLA-[A-Z0-9]+\*\d+:\d+)(?::\d+)*$', '${1}')

    context = ((protein('mhc.a') == protein('current_mhc_a')) &
               (protein('mhc.b') == protein('current_mhc_b')) &
               (pl.col('species') == pl.col('current_species')) &
               (pl.col('mhc.class') == pl.col('current_mhc_class')))
    exact = context & ((pl.col('matches_original') == 'true') | (pl.col('matches_repaired') == 'true'))
    for name, expected in [('matches_context', context), ('exact_tcr_pmhc_match', exact)]:
        if joined.filter(~pl.col(name).is_in(['true', 'false']) |
                         (expected != (pl.col(name) == 'true'))).height:
            raise ValueError(f'Incorrect {name} structure match flag')


def check_archives(root: Path, inventory: pl.DataFrame) -> int:
    """Require every declared archive member and size, and reject unexpected members."""
    if inventory.select('archive', 'member').n_unique() != inventory.height:
        raise ValueError('Duplicate declared archive member')
    count = 0
    for (name,), part in sorted(inventory.group_by('archive'), key=lambda item: item[0]):
        wanted = dict(part.select('member', pl.col('size_bytes').cast(pl.Int64)).iter_rows())
        seen = set()
        with tarfile.open(root / 'data' / name, mode='r|gz') as archive:
            for member in archive:
                if not member.isfile():
                    continue
                if member.name in seen or wanted.get(member.name) != member.size:
                    raise ValueError(f'Unexpected, duplicate or changed asset: {name}/{member.name}')
                seen.add(member.name)
        if seen != set(wanted):
            raise ValueError(f'Missing structure assets in {name}')
        count += len(seen)
    return count


def check_supplemental(root: Path, records: pl.DataFrame, chains: pl.DataFrame,
                       out: Path) -> dict:
    """Check content-addressed supplemental assets and literal current-record proofs."""
    base = root / 'metadata/supplemental'
    assets = _read(base / 'structures.tsv.gz')
    models = _read(base / 'models.tsv.gz')
    links = _read(base / 'record_links.tsv.gz')
    if assets['structure_id'].n_unique() != assets.height:
        raise ValueError('Duplicate supplemental coordinate identity')
    checked = set()
    for (name,), part in sorted(assets.group_by('archive_path'), key=lambda item: item[0]):
        wanted = dict(part.select('archive_member', 'structure_id').iter_rows())
        with tarfile.open(root / name, 'r|gz') as archive:
            for member in archive:
                if not member.isfile() or member.name not in wanted:
                    continue
                stream = archive.extractfile(member)
                if stream is None or hashlib.file_digest(stream, 'sha256').hexdigest() != wanted[member.name]:
                    raise ValueError('Supplemental coordinate content does not match its identity')
                checked.add(wanted[member.name])
    if checked != set(assets['structure_id']):
        raise ValueError('Missing supplemental coordinate assets')
    key = ['structure_id', 'record_id', 'sequence_state']
    if links.select(key).n_unique() != links.height:
        raise ValueError('Duplicate supplemental sequence-state link')
    contexts = records.select('record_id', *[
        pl.col(name).alias('current_' + name) for name in
        ['species', 'antigen.epitope', 'mhc.a', 'mhc.b', 'mhc.class']])
    states = []
    for state, column in [('original', 'cdr3.original'), ('repaired', 'cdr3')]:
        frame = contexts
        for gene, side in [('TRA', 'alpha'), ('TRB', 'beta')]:
            frame = frame.join(chains.filter(pl.col('gene') == gene).select(
                'record_id', pl.col(column).alias('current_' + side)),
                on='record_id', how='left', validate='1:1')
        states.append(frame.with_columns(pl.lit(state).alias('sequence_state')))
    proof = links.join(pl.concat(states), on=['record_id', 'sequence_state'],
                       how='left', validate='m:1').join(models, on='structure_id', validate='m:m')
    def protein(name):
        return pl.col(name).str.replace(r'^H-2-?', 'H2-').str.replace(
            r'^(HLA-[A-Z0-9]+\*\d+:\d+)(?::\d+)*$', '${1}')

    valid = pl.all_horizontal([
        pl.col('current_alpha').fill_null('') != '', pl.col('current_beta').fill_null('') != '',
        pl.col('coordinate_alpha').str.count_matches(pl.col('current_alpha'), literal=True) == 1,
        pl.col('coordinate_beta').str.count_matches(pl.col('current_beta'), literal=True) == 1,
        pl.col('coordinate_peptide') == pl.col('current_antigen.epitope'),
        *[pl.col(name) == pl.col('current_' + name) for name in
          ['species', 'antigen.epitope', 'mhc.class']],
        *[protein(name) == protein('current_' + name) for name in ['mhc.a', 'mhc.b']],
    ])
    proved = proof.filter(valid).select(key).unique()
    if links.join(proved, on=key, how='anti').height:
        raise ValueError('Supplemental record link lacks literal sequence/pMHC proof')
    files = links.join(assets, on='structure_id', how='left', validate='m:1').sort(key)
    if files['archive_path'].null_count():
        raise ValueError('Supplemental linked coordinate is absent')
    files.write_parquet(out / 'supplemental_structure_files.parquet')
    return {'supplemental_structures': assets.height, 'supplemental_links': links.height,
            'supplemental_linked_records': links['record_id'].n_unique()}


def run(root: Path, tables: Path, out: Path, manifest: Path, *, fetch: bool = False) -> dict:
    """Verify the independent mapping, raw files and all declared native accessions."""
    inputs = verify_inputs(root, manifest, fetch=fetch)
    records, chains = read_table(tables, 'records'), read_table(tables, 'chains')
    links = _read(root / 'metadata/current/record_links.tsv.gz')
    source = _read(root / 'vdjdb_structures_metadata.tsv.gz').with_row_index('source_row', offset=1)
    source = source.with_columns(pl.col('source_row').cast(pl.String))
    check_links(records, chains, links, source)
    inventory = _read(root / 'metadata/current/asset_inventory.tsv.gz')
    assets = check_archives(root, inventory)
    mapping = _read(root / 'metadata/current/asset_sources.tsv.gz')
    if (mapping.join(inventory.select('asset_row'), on='asset_row', how='anti').height or
            mapping.join(source.select('source_row'), on='source_row', how='anti').height or
            inventory.join(mapping.select('asset_row').unique(), on='asset_row', how='anti').height):
        raise ValueError('Incomplete or invalid structure asset/source mapping')
    deposited = _read(root / 'data/pdb_sources/manifest.tsv.gz')
    for item in deposited.iter_rows(named=True):
        with gzip.open(root / item['file'], 'rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != item['source_sha256']:
                raise ValueError('Deposited coordinate source checksum mismatch')
    # Publication deduplication can retain a PDB observation under another chunk.
    # Check accessions across all records rather than filtering the aggregate chunk.
    declared = records.filter(pl.col('meta.structure.id').str.contains(r'^[0-9][A-Za-z0-9]{3}$')).select(
        'record_id', pl.col('meta.structure.id').str.to_lowercase().alias('pdb_id'))
    native = inventory.filter(pl.col('archive') == 'pdb_files_native.tgz').select(
        pl.col('asset_pdb_id').str.to_lowercase().alias('pdb_id'), 'archive', 'member')
    coverage = declared.join(native, on='pdb_id', how='left', validate='m:1').join(
        deposited.select('pdb_id', 'file', 'source_sha256'), on='pdb_id', how='left', validate='m:1')
    if coverage.filter(pl.col('member').is_null() & pl.col('file').is_null()).height:
        raise ValueError('Current declared PDB accession lacks deposited coordinates')
    files = inventory.filter(pl.col('archive').is_in(['pdb_files.tgz', 'pdb_files_native.tgz'])).select(
        'asset_row', 'archive', 'member')
    files = files.join(mapping, on='asset_row', validate='1:m').join(links.select(
        'source_row', 'record_id', 'tcr_pmhc_hash', 'link_basis', 'matches_original',
        'matches_repaired', 'exact_tcr_pmhc_match'), on='source_row', validate='m:m').select(
        'record_id', 'source_row', 'tcr_pmhc_hash', 'archive', 'member', 'link_basis',
        'matches_original', 'matches_repaired', 'exact_tcr_pmhc_match').unique().sort(
        'record_id', 'source_row', 'archive', 'member')
    out.mkdir(parents=True, exist_ok=True)
    files.write_parquet(out / 'structure_files.parquet')
    coverage.sort('record_id').write_parquet(out / 'native_structure_files.parquet')
    supplemental = (check_supplemental(root, records, chains, out)
                    if (root / 'metadata/supplemental/structures.tsv.gz').is_file() else {})
    report = {'dataset_revision': REVISION, 'verified_inputs': inputs, 'source_rows': source.height,
              'record_links': links.height, 'linked_records': links['record_id'].n_unique(),
              'archive_assets': assets, 'coordinate_record_links': files.height,
              'declared_native_accessions': declared['pdb_id'].n_unique(),
              'deposited_sources': deposited.height, **supplemental, 'status': 'passed'}
    (out / 'structure_files.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
