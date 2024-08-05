import pandas as pd

COMPLEX_SLIM_ANNOT_COLS = [
    "gene",
    "cdr3",
    "species",
    "antigen.epitope",
    "antigen.gene",
    "antigen.species"
]

SUMMARY_COLS = ["complex.id",
                "v.segm", "j.segm",
                "mhc.a", "mhc.b", "mhc.class",
                "reference.id", "vdjdb.score"]


def generate_slim_db(default_db: pd.DataFrame):
    print("Generating and writing slim database")
    slim_db = default_db[COMPLEX_SLIM_ANNOT_COLS + SUMMARY_COLS]
    slim_db['j.start'] = default_db['cdr3fix'].apply(lambda x: x['jStart'] if not x == '' else x)
    slim_db['v.end'] = default_db['cdr3fix'].apply(lambda x: x['vEnd'] if not x == '' else x)
    slim_db.drop_duplicates(subset=COMPLEX_SLIM_ANNOT_COLS, inplace=True)
    slim_db.to_csv('../database/vdjdb_slim.txt', sep='\t')
