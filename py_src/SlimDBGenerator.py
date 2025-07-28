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
                "reference.id", "vdjdb.score", "vdjdb.pgen.score"]


def _aggregating_function(x):
    return ','.join(sorted(set(x.dropna().apply(str))))


def generate_slim_db(default_db: pd.DataFrame) -> None:
    """
    Slim db generator. Write it to /database/ folder
    :param default_db: default vdj db DataFrame
    """
    slim_db = default_db[COMPLEX_SLIM_ANNOT_COLS + SUMMARY_COLS]
    slim_db['j.start'] = default_db['cdr3fix'].apply(lambda x: x['jStart'] if not x == '' else x)
    slim_db['v.end'] = default_db['cdr3fix'].apply(lambda x: x['vEnd'] if not x == '' else x)
    scores_aggregates = slim_db.groupby(COMPLEX_SLIM_ANNOT_COLS)['vdjdb.score'].max().reset_index()['vdjdb.score']
    slim_db = slim_db.groupby(COMPLEX_SLIM_ANNOT_COLS).agg(_aggregating_function).reset_index()
    slim_db['vdjdb.score'] = scores_aggregates
    slim_db.set_index('gene').to_csv('../database/vdjdb.slim.txt', sep='\t', quotechar='"')
