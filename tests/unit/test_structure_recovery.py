"""Structure association preserves context; optional domain proof checks complete sequences."""
from __future__ import annotations

import io
import tarfile

import polars as pl
import pytest

from vdjdb.assemble.evidence import structures

ALPHA = ('GQNIDQPTEMTATEGAIVQINCTYQTSGFNGLFWYQQHAGEAPTFLSYNVLDGLEEKGRFSSFLSRSKGYSYLLL'
         'KELQMKDSASYLCAVRDSNYQLIWGAGTKLIIKP')
BETA = ('GAVVSQHPSWVICKSGTSVKIECRSLDFQATTMFWYRQFPKQSLMLMATSNEGSKATYEQGVEKDKFLINHASL'
        'TLSTLTVTSAHPEDSSFYICASSLGQETQYFGPGTRLTVT')


def inputs(*, association=False):
    records = pl.DataFrame({'record_id':['r1'], 'species':['HomoSapiens'],
        'antigen.epitope':['PEPTIDE'], 'mhc.a':['HLA-A*02:01'], 'mhc.b':['B2M'],
        'mhc.class':['MHCI']})
    chains = pl.DataFrame({'record_id':['r1','r1'], 'gene':['TRA','TRB'],
        'cdr3':['CAVRDSNYQLIW','CASSLGQETQYF'],
        'v.segm':['TRAV1-2*01','TRBV20-1*01'], 'j.segm':['TRAJ33*01','TRBJ2-7*01'],
        'TCR_hash':['current','current']})
    metadata = records.drop('record_id').with_columns(
        pl.lit('source').alias('tcr_pmhc_hash'), pl.lit('False').alias('is_native'),
        pl.lit('CAVRDSNYQLIW' if association else 'UNMATCHED').alias('cdr3.alpha'),
        pl.lit('CASSLGQETQYF').alias('cdr3.beta'),
        pl.lit('UNKNOWN').alias('v.alpha'), pl.lit('TRAJ33*01').alias('j.alpha'),
        pl.lit('TRBV20-1*01').alias('v.beta'), pl.lit('TRBJ2-7*01').alias('j.beta'),
        pl.lit('0.8').alias('iptm'))
    return records,chains,metadata


def archive(path, *, alpha=ALPHA, beta=BETA, peptide='PEPTIDE', gap=None):
    text = 'Chain\tResidue\tResNum\n'
    for chain, sequence in [('TCR_alpha',alpha),('TCR_beta',beta),('peptide',peptide)]:
        for i,aa in enumerate(sequence,1):
            number=i+int(gap is not None and chain==gap[0] and i>=gap[1])
            text += f'{chain}\t{aa}\t{number}\n'
    body=text.encode()
    with tarfile.open(path,'w:gz') as out:
        entry=tarfile.TarInfo('source_aa_coordinates.tsv')
        entry.size=len(body)
        out.addfile(entry,io.BytesIO(body))
    return path


@pytest.mark.parametrize('boundary',['','E'])
def test_recovery_proves_both_full_domains_and_keeps_source_hash(tmp_path,boundary):
    records,chains,metadata=inputs()
    path=archive(tmp_path/'coordinates.tgz',beta=BETA+boundary)
    assert structures(chains,metadata).is_empty()
    evidence=structures(chains,metadata,records=records,coordinates=path)
    assert evidence['record_id'].to_list()==['r1']
    assert evidence['evidence_source'].to_list()==['source']
    native=structures(chains,metadata.with_columns(pl.lit('True').alias('is_native'),
        pl.lit('1abc').alias('meta.structure.id')),records=records,coordinates=path)
    assert native['evidence_value'].to_list()==['1abc']
    assert 'source' in native['evidence_id'].item()


@pytest.mark.parametrize('change',[
    {'species':'MusMusculus'}, {'antigen.epitope':'OTHER'}, {'mhc.a':'HLA-A*02:02'},
    {'mhc.a':'HLA-A*02:01N'}, {'mhc.b':'OTHER'}, {'mhc.class':'MHCII'},
])
def test_recovery_does_not_change_context(tmp_path,change):
    records,chains,metadata=inputs()
    metadata=metadata.with_columns([pl.lit(v).alias(k) for k,v in change.items()])
    assert structures(chains,metadata,records=records,
        coordinates=archive(tmp_path/'coordinates.tgz')).is_empty()


@pytest.mark.parametrize('change',[{'alpha':ALPHA+'E'},{'beta':BETA[:-1]},
                                   {'beta':BETA+'EE'},{'peptide':'PEPTIDD'},{'beta':''}])
def test_coordinates_do_not_allow_partial_or_incompatible_domains(tmp_path,change):
    records,chains,metadata=inputs()
    assert structures(chains,metadata,records=records,
        coordinates=archive(tmp_path/'coordinates.tgz',**change)).is_empty()


def test_domain_proof_requires_explicit_known_alleles(tmp_path):
    records,chains,metadata=inputs()
    path=archive(tmp_path/'coordinates.tgz')
    for v in ['TRAV1-2','TRAV1-2*99','TRAV1-2*01,TRAV1-2*02']:
        changed=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit(v))
            .otherwise(pl.col('v.segm')).alias('v.segm'))
        assert structures(changed,metadata,records=records,coordinates=path).is_empty()


def test_reviewed_spelling_alias_never_resolves_ambiguous_gene(tmp_path):
    from vdjdb.assemble.structures import recover_links
    records,chains,metadata=inputs(association=True)
    exact=pl.DataFrame(schema={'tcr_pmhc_hash':pl.String,'record_id':pl.String})
    path=archive(tmp_path/'empty.tgz',alpha='NN',beta='NN')
    metadata=metadata.with_columns(pl.lit('TRAV1-2*01').alias('v.alpha'),
        pl.lit('TRBJ2.7*01').alias('j.beta'))
    assert recover_links(records,chains,metadata,path,exact)['record_id'].to_list()==['r1']
    metadata=metadata.with_columns(pl.lit('TRBV8').alias('v.beta'))
    assert recover_links(records,chains,metadata,path,exact).is_empty()


def test_exact_hash_cannot_override_species_or_mhc_class():
    records,chains,metadata=inputs()
    metadata=metadata.with_columns(pl.lit('current').alias('tcr_pmhc_hash'))
    for column,value in [('species','MusMusculus'),('mhc.class','MHCII')]:
        with pytest.raises(ValueError,match=r'species|class'):
            structures(chains,metadata.with_columns(pl.lit(value).alias(column)),records=records)


def test_exact_beta_ending_in_e_is_not_unconditionally_trimmed(tmp_path,monkeypatch):
    from vdjdb.assemble import structures as recovery
    records,chains,metadata=inputs()
    monkeypatch.setattr(recovery,'_domains',lambda *_: pl.DataFrame({
        'record_id':['r1'],'alpha_domain':[ALPHA],'beta_domain':[BETA+'E']}))
    evidence=structures(chains,metadata,records=records,
        coordinates=archive(tmp_path/'coordinates.tgz',beta=BETA+'E'))
    assert evidence['record_id'].to_list()==['r1']


def test_synonymous_hla_fields_preserve_expression_suffixes(tmp_path):
    records,chains,metadata=inputs()
    path=archive(tmp_path/'coordinates.tgz')
    for call in ['HLA-A*02:01:01','HLA-A*02:01:01:01']:
        assert structures(chains,metadata.with_columns(pl.lit(call).alias('mhc.a')),
            records=records,coordinates=path)['record_id'].to_list()==['r1']
    assert structures(chains,metadata.with_columns(pl.lit('HLA-A*02:01:01N').alias('mhc.a')),
        records=records,coordinates=path).is_empty()


def test_full_reference_is_explicit_and_supplements_only_missing_alleles(tmp_path,monkeypatch):
    from vdjtools.model import load_bundled

    from vdjdb.annotate.contig import _flanks
    from vdjdb.assemble import structures as recovery
    records,chains,metadata=inputs()
    pre,_suf=_flanks(load_bundled('TRA','olga',organism='human'),'TRA')
    monkeypatch.setattr(recovery,'full_reference_root',lambda: tmp_path)
    monkeypatch.setattr(recovery,'_full_flanks',lambda *_: (
        {'TRAV1-2*99':pre['TRAV1-2*01'],'TRAV1-2*01':'GAA'},{}))
    path=archive(tmp_path/'coordinates.tgz')
    changed=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit('TRAV1-2*99'))
        .otherwise(pl.col('v.segm')).alias('v.segm'))
    assert structures(changed,metadata,records=records,coordinates=path).is_empty()
    for part in [chains,changed]:
        assert structures(part,metadata,records=records,coordinates=path,
            full_germline=True)['record_id'].to_list()==['r1']


def test_full_reference_rejects_missing_environment_and_raw_inputs(tmp_path,monkeypatch):
    from vdjdb.assemble.structures import _full_flanks
    monkeypatch.delenv('ARDA_HOME',raising=False)
    with pytest.raises(ValueError,match='ARDA_HOME'):
        _full_flanks('TRA','human')
    monkeypatch.setenv('ARDA_HOME',str(tmp_path))
    with pytest.raises(ValueError,match=r'local.*reference'):
        _full_flanks('TRA','human')


def test_cli_full_reference_rejects_missing_environment_before_build(tmp_path,monkeypatch):
    from typer.testing import CliRunner

    from vdjdb.cli import app

    monkeypatch.delenv('ARDA_HOME',raising=False)
    metadata=tmp_path/'input.tsv'
    metadata.write_text('tcr_pmhc_hash\tis_native\n')
    result=CliRunner().invoke(app,['build','--structure-metadata',str(metadata),
        '--structure-coordinates',str(archive(tmp_path/'coordinates.tgz')),
        '--structure-full-germline'])
    assert result.exit_code==2
    assert 'ARDA_HOME' in result.output


@pytest.mark.parametrize('gap',[('TCR_alpha',2),('TCR_beta',2),('peptide',2),('TCR_alpha',1)])
def test_gapped_or_offset_residues_cannot_collapse_to_a_complete_domain(tmp_path,gap):
    records,chains,metadata=inputs()
    # The residue letters still spell the exact target; only their numbering reveals missing sequence.
    evidence=structures(chains,metadata,records=records,
        coordinates=archive(tmp_path/'gapped.tgz',gap=gap))
    assert evidence.is_empty()


def test_pair_association_ignores_source_vj_annotations_and_keeps_source_hash():
    records,chains,metadata=inputs(association=True)
    metadata=metadata.with_columns(pl.lit('TRBV7-9*06').alias('v.beta'),pl.lit('').alias('j.alpha'))
    result=structures(chains,metadata,records=records)
    assert result['record_id'].to_list()==['r1']
    assert result['evidence_source'].to_list()==['source']


@pytest.mark.parametrize('call',['','TRAVNOTKNOWN*01','TRAV1-2*01,TRAV1-2*02',
    'TRAV1-2*01;TRAV2*01'])
def test_pair_association_rejects_missing_and_multicall_current_v(call):
    records,chains,metadata=inputs(association=True)
    chains=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit(call)).otherwise(pl.col('v.segm')).alias('v.segm'))
    assert structures(chains,metadata,records=records).is_empty()


@pytest.mark.parametrize('second_v,expected',[('TRAV1-2*02',1),('',1),('TRAV2*01',0),('TRAVNOTKNOWN*01',0)])
def test_pair_association_counts_v_genes_before_context_filter(second_v,expected):
    records,chains,metadata=inputs(association=True)
    records=pl.concat([records,records.with_columns(pl.lit('r2').alias('record_id'),pl.lit('HLA-A*02:02').alias('mhc.a'))])
    other=chains.with_columns(pl.lit('r2').alias('record_id'),pl.lit('other').alias('TCR_hash'),
        pl.when(pl.col('gene')=='TRA').then(pl.lit(second_v)).otherwise(pl.col('v.segm')).alias('v.segm'))
    assert structures(pl.concat([chains,other]),metadata,records=records).height==expected


@pytest.mark.parametrize('field,value',[('mhc.a','HLA-A*02:02'),('mhc.a','HLA-A*02:01N'),
    ('mhc.b','OTHER'),('species','MusMusculus'),('mhc.class','MHCII')])
def test_pair_association_preserves_provided_context(field,value):
    records,chains,metadata=inputs(association=True)
    assert structures(chains,metadata.with_columns(pl.lit(value).alias(field)),records=records).is_empty()


def test_native_association_provenance_is_pdb_specific_and_rejects_context_conflicts():
    records,chains,metadata=inputs(association=True)
    records=records.with_columns(pl.lit('1abc').alias('meta.structure.id'))
    metadata=metadata.with_columns(pl.lit('true').alias('is_native'),pl.lit('HLA-A*02:792N').alias('mhc.a'),
        pl.lit('1abc').alias('meta.structure.id'))
    unrelated=metadata.with_columns(pl.lit('2abc').alias('meta.structure.id'))
    result=structures(chains,pl.concat([metadata,unrelated]),records=records)
    assert result['evidence_value'].to_list()==['1abc']
    assert 'source' in result['evidence_id'].item()
    records=pl.concat([records,records.with_columns(pl.lit('r2').alias('record_id'),pl.lit('HLA-A*02:02').alias('mhc.a'))])
    chains=pl.concat([chains,chains.with_columns(pl.lit('r2').alias('record_id'),pl.lit('other').alias('TCR_hash'))])
    assert structures(chains,metadata,records=records).is_empty()


def test_native_association_compares_documented_h2_prefix_and_hla_synonymous_fields():
    records,chains,metadata=inputs(association=True)
    metadata=metadata.with_columns(pl.lit('true').alias('is_native'),pl.lit('1abc').alias('meta.structure.id'),
        pl.lit('HLA-A*02:01:01').alias('mhc.a'))
    assert structures(chains,metadata,records=records)['evidence_value'].to_list()==['1abc']
    records=records.with_columns(pl.lit('MusMusculus').alias('species'),pl.lit('H2-Db').alias('mhc.a'))
    chains=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit('TRAV1*01'))
        .otherwise(pl.lit('TRBV13-2*01')).alias('v.segm'))
    metadata=metadata.with_columns(pl.lit('MusMusculus').alias('species'),pl.lit('H-2Db').alias('mhc.a'))
    assert structures(chains,metadata,records=records)['evidence_value'].to_list()==['1abc']


def test_pair_association_accepts_one_shared_delta_v_label_on_current_alpha():
    records,chains,metadata=inputs(association=True)
    chains=chains.with_columns(pl.when(pl.col('gene')=='TRA').then(pl.lit('TRDV1*01')).otherwise(pl.col('v.segm')).alias('v.segm'))
    assert structures(chains,metadata,records=records)['record_id'].to_list()==['r1']


def test_pair_association_preserves_missing_current_v_when_other_candidate_has_unique_gene():
    records,chains,metadata=inputs(association=True)
    records=pl.concat([records,records.with_columns(pl.lit('r2').alias('record_id'))])
    other=chains.with_columns(pl.lit('r2').alias('record_id'),pl.lit('other').alias('TCR_hash'),
        pl.lit('').alias('v.segm'))
    combined=pl.concat([chains,other])
    assert structures(combined,metadata,records=records)['record_id'].to_list()==['r1','r2']
    assert other['v.segm'].to_list()==['','']


def test_native_pdb_provenance_does_not_require_source_mhc_columns():
    records,chains,metadata=inputs(association=True)
    records=records.with_columns(pl.lit('1abc').alias('meta.structure.id'))
    metadata=metadata.with_columns(pl.lit('true').alias('is_native'),pl.lit('1abc').alias('meta.structure.id'))
    metadata=metadata.select('tcr_pmhc_hash','is_native','cdr3.alpha','cdr3.beta','antigen.epitope','meta.structure.id')
    assert structures(chains,metadata,records=records)['evidence_value'].to_list()==['1abc']


def test_pair_association_requires_species_specific_v_authority():
    records,chains,metadata=inputs(association=True)
    records=records.with_columns(pl.lit('UnknownSpecies').alias('species'))
    metadata=metadata.with_columns(pl.lit('UnknownSpecies').alias('species'))
    assert structures(chains,metadata,records=records).is_empty()
