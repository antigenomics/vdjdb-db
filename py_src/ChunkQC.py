import re
import pandas as pd
from collections import defaultdict


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

speciesList = ["homosapiens", "musmusculus", "rattusnorvegicus", "macacamulatta"]



def is_aa_seq_valid(aa_seq: str) -> bool:
    if pd.isnull(aa_seq):
        return True
    return len(aa_seq) > 3 and bool(re.match(r'^[ARNDCQEGHILKMFPSTWYV]+$', aa_seq))

def is_MHC_valid(hla_allele: str) -> bool:
    return bool(re.match(r'^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$', hla_allele)) or hla_allele[0:3] != 'HLA'

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
    'reference.id': lambda x: x.startswith('PMID:') or x.startswith('doi:') or x.startswith('http://')
                              or x.startswith('https://') or 'unpublished' in x.lower() if not pd.isnull(x) else True

}


class ChunkQC:
    def __init__(self, chunk_df: pd.DataFrame):

        self.chunk_df = chunk_df

    def check_exist(self):
        if not len(self.chunk_df):
            raise ValueError('Empty file')

    def check_header(self):
        if len(self.chunk_df.columns) != len(set(self.chunk_df.columns)):
            raise ValueError(f'Duplicate columns found: {self.chunk_df.columns}')

        missing_columns = set(ALL_COLS).difference(set(self.chunk_df.columns))
        if missing_columns:
            raise ValueError(f'The following columns are missing: {missing_columns}')

    def process_chunk(self) -> dict[int: str]:

        chunk_error_messages = defaultdict(list)

        self.check_exist()
        self.check_header()

        chunk_duplicates = self.chunk_df[self.chunk_df.duplicated(subset=SIGNATURE_COLS)]
        if len(chunk_duplicates):
            for duplicated_row_ind in chunk_duplicates.index:
                chunk_error_messages[duplicated_row_ind].append('duplicate')

        for validating_column, validator in zip(validators.keys(), validators.values()):
            validator_res_mask = self.chunk_df[validating_column].apply(lambda x: validator(x))
            if not all(validator_res_mask):
                for broken_row_ind in self.chunk_df[~validator_res_mask][SIGNATURE_COLS].index:
                    chunk_error_messages[broken_row_ind].append(f'bad {validating_column}')

        empty_cdr3_rows = self.chunk_df[self.chunk_df.T.apply(lambda x: pd.isnull(x["cdr3.alpha"])
                                                              and pd.isnull(x["cdr3.beta"]))][SIGNATURE_COLS]
        for empty_row_ind in empty_cdr3_rows.index:
            chunk_error_messages[tuple(empty_row_ind)].append('no.cdr3')

        empty_epitope_rows = self.chunk_df[self.chunk_df["antigen.epitope"].apply(pd.isnull)][SIGNATURE_COLS]
        for empty_row_ind in empty_epitope_rows.index:
            chunk_error_messages[tuple(empty_row_ind)].append('no.antigen.seq')

        empty_mhc_rows = self.chunk_df[self.chunk_df.T.apply(lambda x: pd.isnull(x["mhc.a"])
                                                             or pd.isnull(x["mhc.b"]))][SIGNATURE_COLS]
        for empty_row_index in empty_mhc_rows.index:
            chunk_error_messages[tuple(empty_row_index)].append('no.mhc')

        return chunk_error_messages


