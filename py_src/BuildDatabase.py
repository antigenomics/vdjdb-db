import pandas as pd
import argparse
import os
import re
import warnings

from collections import defaultdict

# Database specification

COMPLEX_COLUMNS = [
    "cdr3.alpha",
    "v.alpha",
    "j.alpha",
    "cdr3.beta",
    "v.beta",
    "d.beta",
    "j.beta",
    "species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
    "antigen.epitope",
    "antigen.gene",
    "antigen.species",
    "reference.id"
]

METHOD_COLUMNS = [
    "method.identification",
    "method.frequency",
    "method.singlecell",
    "method.sequencing",
    "method.verification"
]

META_COLUMNS = [
    "meta.study.id",
    "meta.cell.subset",
    "meta.subject.cohort",
    "meta.subject.id",
    "meta.replica.id",
    "meta.clone.id",
    "meta.epitope.id",
    "meta.tissue",
    "meta.donor.MHC",
    "meta.donor.MHC.method",
    "meta.structure.id"
]

ALL_COLS = COMPLEX_COLUMNS + METHOD_COLUMNS + META_COLUMNS

SIGNATURE_COLS = [
    "cdr3.alpha",
    "v.alpha",
    "j.alpha",
    "cdr3.beta",
    "v.beta",
    "d.beta",
    "j.beta",
    "species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
    "antigen.epitope",
    "antigen.gene",
    "antigen.species",
    "reference.id",
    "meta.study.id",
    "meta.cell.subset",
    "meta.subject.cohort",
    "meta.subject.id",
    "meta.replica.id",
    "meta.clone.id",
    "meta.tissue"
]

#  Misc utils and classes, temporary files

os.makedirs('../tmp/', exist_ok=True)


# Table utils

def check_header(header: list):
    if not header:
        raise ValueError('Empty file')
    if len(header) != len(set(header)):
        raise ValueError(f'Duplicate columns found: {header}')

    missing_columns = set(ALL_COLS).difference(set(header))
    if missing_columns:
        raise ValueError(f'The following columns are missing: {missing_columns}')


def read_chunk(path) -> pd.DataFrame:
    chunk_df = pd.read_csv(path, sep='\t', encoding_errors='ignore')
    check_header(list(chunk_df.columns))
    return chunk_df


# Basic validation utils

def is_aa_seq_valid(aa_seq: str) -> bool:
    if pd.isnull(aa_seq):
        return True
    return len(aa_seq) > 3 and bool(re.match(r'^[ARNDCQEGHILKMFPSTWYV]+$', aa_seq))


def is_MHC_valid(hla_allele: str) -> bool:
    return bool(re.match(r'^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$', hla_allele)) or hla_allele[0:3] != 'HLA'


speciesList = ["homosapiens", "musmusculus", "rattusnorvegicus", "macacamulatta"]

validators = {
    "cdr3.alpha": is_aa_seq_valid,
    "v.alpha": lambda x: x.startswith("TRAV") if not pd.isnull(x) else True,
    'j.alpha': lambda x: x.startswith('TRAJ') if not pd.isnull(x) else True,
    'cdr3.beta': is_aa_seq_valid,
    "v.beta": lambda x: x.startswith('TRBV') if not pd.isnull(x) else True,
    'd.beta': lambda x: x.startswith('TRBD') if not pd.isnull(x) else True,
    'j.beta': lambda x: x.startswith('TRBJ') if not pd.isnull(x) else True,
    'species': lambda x: x.lower() in speciesList,
    'mhc.a': is_MHC_valid,
    'mhc.b': is_MHC_valid,
    'mhc.class': lambda x: x == 'MHCI' or x == 'MHCII',
    'antigen.epitope': is_aa_seq_valid,
    'reference.id': lambda x: x.startswith('PMID:') or x.startswith('doi:') or x.startswith('http://"')
                              or x.startswith('https://') or 'unpublished' in x.lower() if not pd.isnull(x) else True

}

antigen_df = pd.read_csv("../patches/antigen_epitope_species_gene.dict", sep='\t', index_col=0)
aggregated_species = antigen_df.groupby(level=0)['antigen.species'].agg(lambda x: x.iloc[0] if len(x) == 1 else list(x))
aggregated_gene = antigen_df.groupby(level=0)['antigen.gene'].agg(lambda x: x.iloc[0] if len(x) == 1 else list(x))

# Chunks reading

parser = argparse.ArgumentParser(description='Arguments for database building')
parser.add_argument('--chunks-to-build', help='chunks to include in database', nargs='+', type=str)
parser.add_argument('--no2fix', action='store_true', help='Fix did not occurred if enabled')
args = parser.parse_args()

chunk_files = set(os.listdir('../chunks')).intersection(args.chunks_to_build) if \
    args.chunks_to_build else os.listdir('../chunks')

chunk_files = [chunk_file for chunk_file in chunk_files if chunk_files[0] != '.']

if not chunk_files:
    raise FileNotFoundError('No database chunks to process')

chunk_df_list = []

for chunk_file in chunk_files:

    print(f'processing chunk {chunk_file}')
    chunk_df = read_chunk(f'../chunks/{chunk_file}')
    chunk_error_messages = defaultdict(list)
    chunk_duplicates = chunk_df[chunk_df.duplicated(subset=SIGNATURE_COLS)]
    if len(chunk_duplicates):
        for duplicated_row_ind in chunk_duplicates.index:
            chunk_error_messages[duplicated_row_ind].append('duplicate')

    for validating_column, validator in zip(validators.keys(), validators.values()):
        validator_res_mask = chunk_df[validating_column].apply(lambda x: validator(x))
        if not all(validator_res_mask):
            for broken_row_ind in chunk_df[~validator_res_mask][SIGNATURE_COLS].index:
                chunk_error_messages[broken_row_ind].append(f'bad {validating_column}')

    empty_cdr3_rows = chunk_df[chunk_df.T.apply(lambda x: pd.isnull(x["cdr3.alpha"])
                                                          and pd.isnull(x["cdr3.beta"]))][SIGNATURE_COLS]
    for empty_row_ind in empty_cdr3_rows.index:
        chunk_error_messages[tuple(empty_row_ind)].append('no.cdr3')

    empty_epitope_rows = chunk_df[chunk_df["antigen.epitope"].apply(pd.isnull)][SIGNATURE_COLS]

    for empty_row_ind in empty_epitope_rows.index:
        chunk_error_messages[tuple(empty_row_ind)].append('no.antigen.seq')

    empty_mhc_rows = chunk_df[chunk_df.T.apply(lambda x: pd.isnull(x["mhc.a"])
                                                         or pd.isnull(x["mhc.b"]))][SIGNATURE_COLS]

    for empty_row_index in empty_mhc_rows.index:
        chunk_error_messages[tuple(empty_row_index)].append('no.mhc')

    if chunk_error_messages.keys():
        print(dict(chunk_error_messages))
        warnings.warn(f"There were errors processing {chunk_file}")
        continue

    chunk_df['antigen.species'] = chunk_df['antigen.epitope'].apply(
    lambda x: aggregated_species[x] if x in aggregated_species else None
)
    chunk_df['antigen.gene'] = chunk_df['antigen.epitope'].apply(
    lambda x: aggregated_gene[x] if x in aggregated_gene else None
)
    chunk_df_list.append(chunk_df)

# Fix CDR3 sequences
# this part is skipped

os.makedirs('../database/', exist_ok=True)
pd.concat(chunk_df_list).to_csv('../database/vdjdb.txt', sep='\t')
