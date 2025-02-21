import pandas as pd
import json
import csv
from ChunkQC import SIGNATURE_COLS, METHOD_COLUMNS, META_COLUMNS
import sys
from multiprocessing import Pool
import numpy as np

sys.path.append('../../')
sys.path.append('../../mirpy')

from mirpy.mir.basic import pgen

olga_pgen_human_trb = pgen.OlgaModel(model='../../mirpy/mir/resources/olga/default_models/human_T_beta')
olga_pgen_human_tra = pgen.OlgaModel(model='../../mirpy/mir/resources/olga/default_models/human_T_alpha',
                                     is_d_present=False)
olga_pgen_mouse_trb = pgen.OlgaModel(model='../../mirpy/mir/resources/olga/default_models/mouse_T_beta')
olga_pgen_mouse_tra = pgen.OlgaModel(model='../../mirpy/mir/resources/olga/default_models/mouse_T_alpha',
                                     is_d_present=False)


models_dict = {
    'homosapiens_TRB': olga_pgen_human_trb,
    'homosapiens_TRA': olga_pgen_human_tra,
    'musmusculus_TRB': olga_pgen_mouse_trb,
    'musmusculus_TRA': olga_pgen_mouse_tra
}

VERY_HIGH_CONFIDENCE_CUTOFF_B = -7.3
HIGH_CONFIDENCE_CUTOFF_B = -12.1
MEDIUM_CONFIDENCE_CUTOFF_B = -15.6

VERY_HIGH_CONFIDENCE_CUTOFF_A = -5.7
HIGH_CONFIDENCE_CUTOFF_A = -10.1
MEDIUM_CONFIDENCE_CUTOFF_A = -15.1

VERY_HIGH_CONFIDENCE_CUTOFF_MOUSE_B = -5.4
HIGH_CONFIDENCE_CUTOFF_MOUSE_B = -6.6
MEDIUM_CONFIDENCE_CUTOFF_MOUSE_B = -8.2

VERY_HIGH_CONFIDENCE_CUTOFF_MOUSE_A = -4.5
HIGH_CONFIDENCE_CUTOFF_MOUSE_A = -8.3
MEDIUM_CONFIDENCE_CUTOFF_MOUSE_A = -10

def calc_pgen(multiargument):
    cdr3aa = multiargument[0]
    gene = multiargument[1]
    specie = multiargument[2].lower()
    if specie not in {'musmusculus', 'homosapiens'}:
        return None
    model = models_dict[f'{specie}_{gene}']
    p_gen = model.compute_pgen_cdr3aa(cdr3aa)
    log10_pgen = np.log10(p_gen)
    return log10_pgen


def get_web_method(method_identification: str) -> str:
    method_identification = method_identification.lower()
    if "sort" in method_identification:
        return "sort"
    elif "culture" in method_identification or "cloning" in method_identification or "targets" in method_identification:
        return "culture"
    else:
        return "other"


def get_web_method_seq(clone_row: pd.Series) -> str:
    if clone_row["method.singlecell"] != "":
        return "singlecell"
    else:
        method_data = clone_row["method.sequencing"].lower()
        if "sanger" in method_data:
            return "sanger"
        elif "-seq" in method_data:
            return "amplicon"
        else:
            return "other"


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


def generate_default_db(master_table: pd.DataFrame) -> pd.DataFrame:
    """
    Generates vdjdb default txt file from full table. Writes it to /database/ folder
    :param master_table: full vdj db table
    :return: default vdj db table
    """
    master_table.fillna("", inplace=True)
    complex_id_count = 0

    sample_counts = master_table.value_counts(subset=SIGNATURE_COLS_PER_SAMPLE)
    study_counts = master_table.set_index(SIGNATURE_COLS)
    clones_list = []

    for _, clone in master_table.iterrows():

        if not (clone["cdr3.alpha"] == "") or (clone["cdr3.beta"] == ""):
            complex_id_count += 1
            complex_id = complex_id_count
        else:
            complex_id = 0

        for chain in ["alpha", "beta"]:
            if clone[f"cdr3.{chain}"] and clone[f"v.{chain}"] and clone[f"j.{chain}"]:
                clone_compact = {
                    "complex.id": complex_id,
                    "gene": "TRA" if chain == "alpha" else "TRB",
                    "cdr3": clone[f"cdr3.{chain}"],
                    "v.segm": clone[f"v.{chain}"],
                    "j.segm": clone[f"j.{chain}"],
                }

                for coll in COMPLEX_ANNOT_COLS:
                    clone_compact[coll] = clone[coll]

                clone_compact["method"] = {coll.split("method.")[1]: clone[coll] for coll in METHOD_COLUMNS}
                clone_compact["meta"] = {coll.split("meta.")[1]: clone[coll] for coll in META_COLUMNS}
                clone_compact["meta"]["samples.found"] = int(sample_counts.loc[tuple(clone[coll] for coll
                                                                                 in SIGNATURE_COLS_PER_SAMPLE)])
                clone_compact["meta"]["studies.found"] \
                    = len(study_counts.loc[tuple(clone[coll]
                                                 for coll in
                                                 SIGNATURE_COLS)].index.get_level_values(14).unique())

                clone_compact["cdr3fix"] = clone[f"cdr3fix.{chain}"]
                clone_compact["web.method"] = get_web_method(clone["method.identification"])
                clone_compact["web.method.seq"] = get_web_method_seq(clone)
                if clone_compact["cdr3fix"] != "":
                    clone_compact["web.cdr3fix.nc"] = "no" if clone_compact["cdr3fix"]["jCanonical"] \
                                                              and clone_compact["cdr3fix"]["vCanonical"] else "yes"
                    clone_compact["web.cdr3fix.unmp"] = "no" if clone_compact["cdr3fix"]["vEnd"] > -1 \
                                                                and clone_compact["cdr3fix"]["jStart"] else "yes"

                clones_list.append(clone_compact)
    default_db = pd.DataFrame(clones_list)

    with Pool(24) as p:
        log_10_pgens = list(p.map(calc_pgen, [(x.cdr3, x.gene, x.species) for _, x in default_db.iterrows()]))

    default_db['log_10_pgen'] = log_10_pgens

    homosapiens_beta_score = pd.cut(default_db[(default_db.species == 'HomoSapiens') & (default_db.gene == 'TRB')].log_10_pgen,
                              [-500, MEDIUM_CONFIDENCE_CUTOFF_B, HIGH_CONFIDENCE_CUTOFF_B,
                               VERY_HIGH_CONFIDENCE_CUTOFF_B, 0], labels=[0, 1, 2, 3], )
    homosapiens_alpha_score = pd.cut(default_db[(default_db.species == 'HomoSapiens') & (default_db.gene == 'TRA')].log_10_pgen,
                               [-500, MEDIUM_CONFIDENCE_CUTOFF_A, HIGH_CONFIDENCE_CUTOFF_A,
                                VERY_HIGH_CONFIDENCE_CUTOFF_A, 0], labels=[0, 1, 2, 3], )
    musmusculus_beta_score = pd.cut(default_db[(default_db.species == 'MusMusculus') & (default_db.gene == 'TRB')].log_10_pgen,
                              [-500, MEDIUM_CONFIDENCE_CUTOFF_MOUSE_B, HIGH_CONFIDENCE_CUTOFF_MOUSE_B,
                               VERY_HIGH_CONFIDENCE_CUTOFF_MOUSE_B, 0], labels=[0, 1, 2, 3], )
    musmusculus_alpha_score = pd.cut(default_db[(default_db.species == 'MusMusculus') & (default_db.gene == 'TRA')].log_10_pgen,
                               [-500, MEDIUM_CONFIDENCE_CUTOFF_A, HIGH_CONFIDENCE_CUTOFF_A,
                                VERY_HIGH_CONFIDENCE_CUTOFF_A, 0], labels=[0, 1, 2, 3], )
    vdj_db_score = pd.concat([homosapiens_beta_score, homosapiens_alpha_score,
                              musmusculus_beta_score, musmusculus_alpha_score])
    vdj_db_score.name = 'vdjdb.score'
    default_db['vdjdb.score'] = vdj_db_score
    default_db['vdjdb.score'] = default_db['vdjdb.score'].fillna(0)
    default_db['vdjdb.score'] = default_db['vdjdb.score'].apply(int)
    default_db.to_pickle('../database/vdjvd.pkl')
    default_db = default_db.drop('log_10_pgen', axis=1) # delete after front fix
    default_db_to_write = default_db.copy()

    for complex_col in ["method", "meta", "cdr3fix"]:
        default_db_to_write[complex_col] = default_db_to_write[complex_col].apply(lambda x: json.dumps(x))

    default_db_to_write.set_index("complex.id").to_csv("../database/vdjdb.txt", sep="\t", quoting=csv.QUOTE_NONE)
    return default_db
