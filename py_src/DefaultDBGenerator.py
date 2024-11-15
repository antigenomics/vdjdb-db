import pandas as pd
import json
import csv
from ChunkQC import SIGNATURE_COLS, METHOD_COLUMNS, META_COLUMNS
import sys
from multiprocessing import Pool

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


def calc_pgen(multiargument):
    cdr3aa = multiargument[0]
    gene = multiargument[1]
    specie = multiargument[2].lower()
    if specie not in {'musmusculus', 'homosapiens'}:
        return None
    model = models_dict[f'{specie}_{gene}']
    p_gen = model.compute_pgen_cdr3aa(cdr3aa)
    return p_gen


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
            if clone[f"cdr3.{chain}"]:
                clone_compact = {
                    "complex.id": complex_id,
                    "gene": "TRA" if chain == "alpha" else "TRB",
                    "cdr3": clone[f"cdr3.{chain}"],
                    "v.segm": clone[f"v.{chain}"],
                    "j.segm": clone[f"j.{chain}"],
                #    "pgen": calc_pgen(clone[f"cdr3.{chain}"], chain, clone['species'].lower())
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
                clone_compact["vdjdb.score"] = 1  #placeholder for tests
                clone_compact["web.method"] = get_web_method(clone["method.identification"])
                clone_compact["web.method.seq"] = get_web_method_seq(clone)
                if clone_compact["cdr3fix"] != "":
                    clone_compact["web.cdr3fix.nc"] = "no" if clone_compact["cdr3fix"]["jCanonical"] \
                                                              and clone_compact["cdr3fix"]["vCanonical"] else "yes"
                    clone_compact["web.cdr3fix.unmp"] = "no" if clone_compact["cdr3fix"]["vEnd"] > -1 \
                                                                and clone_compact["cdr3fix"]["jStart"] else "yes"

                clones_list.append(clone_compact)
    default_db = pd.DataFrame(clones_list)
    for complex_col in ["method", "meta", "cdr3fix"]:
        default_db[complex_col] = default_db[complex_col].apply(lambda x: json.dumps(x))

    with Pool(24) as p:
        pgens = list(p.map(calc_pgen, [(x.cdr3, x.gene, x.species) for _, x in default_db.iterrows()]))

    default_db['pgen'] = pgens

    default_db.set_index("complex.id").to_csv("../database/vdjdb.txt", sep="\t", quoting=csv.QUOTE_NONE)
    return pd.DataFrame(clones_list)
