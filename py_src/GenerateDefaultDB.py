import pandas as pd
from ChunkQC import SIGNATURE_COLS, METHOD_COLUMNS, META_COLUMNS


def get_web_method(method_identification: str) -> str:
    method_identification = method_identification.lower()
    if 'sort' in method_identification:
        return "sort"
    elif 'culture' in method_identification or 'cloning' in method_identification or 'targets' in method_identification:
        return 'culture'
    else:
        return "other"


def get_web_method_seq(clone_row: pd.Series) -> str:
    if not pd.isnull(clone_row["method.singlecell"]):
        return "singlecell"
    else:
        method_data = clone_row["method.sequencing"].lower()
        if 'sanger' in method_data:
            return 'sanger'
        elif '-seq' in method_data:
            return 'amplicon'
        else:
            return 'other'


COMPLEX_ANNOT_COLS = [
    "species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
    "antigen.epitope",
    "antigen.gene",
    "antigen.species",
    "reference.id"]

SIGNATURE_COLS_PER_SAMPLE = [
    "cdr3.alpha",
    "v.alpha",
    "j.alpha",
    "cdr3.beta",
    "v.beta",
    "j.beta",
    "species",
    "mhc.a",
    "mhc.b",
    "mhc.class",
    "antigen.epitope"
]


def generate_default_db(master_table: pd.DataFrame):
    complex_id_count = 0

    sample_counts = master_table.value_counts(subset=SIGNATURE_COLS_PER_SAMPLE)
    study_counts = master_table.set_index(SIGNATURE_COLS)
    clones_list = []

    for _, clone in master_table.iterrows():

        if not (pd.isnull(clone["cdr3.alpha"]) or pd.isnull(clone["cdr3.beta"])):
            complex_id_count += 1
            complex_id = complex_id_count
        else:
            complex_id = 0

        for chain in ['alpha', 'beta']:
            if not pd.isnull(clone[f'cdr3.{chain}']):
                clone_compact = {
                    'complex.id': complex_id,
                    'gene': 'TRA' if chain == 'alpha' else 'TRB',
                    'cdr3': clone[f"cdr3.{chain}"],
                    'v.segm': clone[f"v.{chain}"],
                    'j.segm': clone[f"j.{chain}"], }

                for coll in COMPLEX_ANNOT_COLS:
                    clone_compact[coll] = clone[coll]

                clone_compact['method'] = {coll.split('method.')[1]: clone[coll] for coll in METHOD_COLUMNS}
                clone_compact['meta'] = {coll.split('meta.')[1]: clone[coll] for coll in META_COLUMNS}
                clone_compact['meta']['samples.found'] = sample_counts.loc[(clone[coll] for coll
                                                                            in SIGNATURE_COLS_PER_SAMPLE)]
                clone_compact['meta']['studies.found'] = len(study_counts.loc[(clone[coll]
                                                                               for coll in
                                                                               SIGNATURE_COLS_PER_SAMPLE)][
                                                                 'reference.id'].unique())

                clone_compact['cdr3fix'] = ''
                clone_compact['vdjdb.score'] = 1
                clone_compact['web.method'] = get_web_method(clone['method.identification'])
                clone_compact['web.method.seq'] = get_web_method_seq(clone)
                clone_compact['web.cdr3fix.nc'] = 'no'
                clone_compact['web.cdr3fix.unmp'] = 'no'

                clones_list.append(clone_compact)
    pd.concat(clones_list).to_csv('../database/vdjdb.txt', sep='\t')
    return pd.concat(clones_list)
