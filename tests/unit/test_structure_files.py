"""Independent structure-file audit rejects changed and missing raw inputs."""
import gzip
import hashlib
import io
import tarfile

import polars as pl
import pytest

from vdjdb.validate.structure_files import check_archives, check_links, check_supplemental, verify_inputs


def test_streamed_archive_members_and_pinned_inputs(tmp_path):
    (tmp_path / 'data').mkdir()
    path = tmp_path / 'data/coordinates.tgz'
    with tarfile.open(path, 'w:gz') as archive:
        item = tarfile.TarInfo('one.pdb')
        item.size = 4
        archive.addfile(item, io.BytesIO(b'ATOM'))
    inventory = pl.DataFrame({'archive': ['coordinates.tgz'], 'member': ['one.pdb'], 'size_bytes': ['4']})
    assert check_archives(tmp_path, inventory) == 1
    for invalid in [inventory.with_columns(pl.lit('two.pdb').alias('member')),
                    inventory.with_columns(pl.lit('5').alias('size_bytes'))]:
        with pytest.raises(ValueError):
            check_archives(tmp_path, invalid)
    manifest = tmp_path / 'inputs.tsv'
    pl.DataFrame({'path': ['data/coordinates.tgz'],
                  'sha256': [hashlib.sha256(path.read_bytes()).hexdigest()]}).write_csv(
                      manifest, separator='\t')
    assert verify_inputs(tmp_path, manifest) == 1
    path.write_bytes(b'changed')
    with pytest.raises(ValueError, match='checksum'):
        verify_inputs(tmp_path, manifest)
    path.unlink()
    with pytest.raises(ValueError, match='--fetch'):
        verify_inputs(tmp_path, manifest)


def test_links_recompute_source_identity_and_context_flags():
    record = {'record_id': 'r1', 'species': 'HomoSapiens', 'antigen.epitope': 'PEP',
              'mhc.a': 'HLA-A*02:01', 'mhc.b': 'B2M', 'mhc.class': 'MHCI',
              'reference.id': 'P1', 'meta.structure.id': '', 'chunk.file': 'paper.tsv', 'chunk.row': '1'}
    records = pl.DataFrame([record])
    chains = pl.DataFrame({'record_id': ['r1', 'r1'], 'gene': ['TRA', 'TRB'],
                          'cdr3': ['CAVF', 'CASSF'], 'cdr3.original': ['CAVF', 'CASSF']})
    source = {key: record[key] for key in ['species', 'antigen.epitope', 'mhc.a', 'mhc.b', 'mhc.class']}
    source.update({'source_row': '1', 'tcr_pmhc_hash': 'model1', 'is_native': 'False',
                   'meta.structure.id': '', 'cdr3.alpha': 'CAVF', 'cdr3.beta': 'CASSF'})
    link = {'record_id': 'r1', 'source_row': '1', 'tcr_pmhc_hash': 'model1', 'is_native': 'False',
            'meta.structure.id': '', **{'current_' + key.replace('.', '_'): value
                                       for key, value in record.items() if key != 'record_id'},
            'current_repaired_alpha': 'CAVF', 'current_original_alpha': 'CAVF',
            'current_repaired_beta': 'CASSF', 'current_original_beta': 'CASSF',
            'matches_original': 'true', 'matches_repaired': 'true', 'matches_context': 'true',
            'exact_tcr_pmhc_match': 'true'}
    check_links(records, chains, pl.DataFrame([link]), pl.DataFrame([source]))
    for name, value in [('tcr_pmhc_hash', 'broken'), ('is_native', 'True'),
                        ('current_chunk_row', '2'),
                        ('exact_tcr_pmhc_match', 'not-a-boolean'), ('matches_context', 'false'),
                        ('matches_original', 'invalid')]:
        with pytest.raises(ValueError):
            check_links(records, chains, pl.DataFrame([{**link, name: value}]), pl.DataFrame([source]))


@pytest.mark.parametrize("pairing,disposition", [("reported_assignment", "predicted"),
                                               ("constructed_mismatch", "mock"),
                                               ("unknown", "native")])
def test_supplemental_coordinates_and_record_proof(tmp_path, pairing, disposition):
    base = tmp_path / 'metadata/supplemental'
    base.mkdir(parents=True)
    (tmp_path / 'data').mkdir()
    raw = b'ATOM\n'
    digest = hashlib.sha256(raw).hexdigest()
    with tarfile.open(tmp_path / 'data/extra.tar.gz', 'w:gz') as archive:
        member = tarfile.TarInfo(digest + '.pdb')
        member.size = len(raw)
        archive.addfile(member, io.BytesIO(raw))
    assets = pl.DataFrame({'structure_id': [digest], 'archive_path': ['data/extra.tar.gz'],
                          'archive_member': [digest + '.pdb'], 'coordinate_bytes': [len(raw)]})
    links = pl.DataFrame({'structure_id': [digest], 'record_id': ['r1'], 'sequence_state': ['original'],
                         'link_basis': ['literal_sequence_context']})
    record = {'record_id': 'r1', 'species': 'HomoSapiens', 'antigen.epitope': 'PEP',
              'mhc.a': 'HLA-A*02:01', 'mhc.b': 'B2M', 'mhc.class': 'MHCI'}
    models = pl.DataFrame([{'structure_id': digest, 'pairing_status': pairing,
        'source_disposition': disposition, **{k: v for k, v in record.items() if k != 'record_id'},
        'coordinate_alpha': 'AACAVFAA', 'coordinate_beta': 'AACASSFAA', 'coordinate_peptide': 'PEP'}])
    for name, frame in [('structures', assets), ('models', models), ('record_links', links)]:
        with gzip.open(base / (name + '.tsv.gz'), 'wb') as stream:
            frame.write_csv(stream, separator='\t')
    chains = pl.DataFrame({'record_id': ['r1', 'r1'], 'gene': ['TRA', 'TRB'],
                          'cdr3': ['CAVF', 'CASSF'], 'cdr3.original': ['CAVF', 'CASSF']})
    assert check_supplemental(tmp_path, pl.DataFrame([record]), chains, tmp_path)['supplemental_links'] == 1
    with pytest.raises(ValueError, match='literal'):
        check_supplemental(tmp_path, pl.DataFrame([{**record, 'mhc.a': 'HLA-A*01:01'}]), chains, tmp_path)
    for side, junction in [('alpha', 'CAVF'), ('beta', 'CASSF')]:
        with gzip.open(base / 'models.tsv.gz', 'wb') as stream:
            models.with_columns(pl.lit(junction * 2).alias('coordinate_' + side)).write_csv(
                stream, separator='\t')
        with pytest.raises(ValueError, match='literal'):
            check_supplemental(tmp_path, pl.DataFrame([record]), chains, tmp_path)
