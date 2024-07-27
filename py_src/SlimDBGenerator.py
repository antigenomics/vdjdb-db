import pandas as pd

slim_db_columns = ['antigen.epitope',
                   'antigen.gene',
                   'antigen.species',
                   'cdr3',
                   'complex.id',
                   'gene',
                   'j.segm',
                   'mhc.a',
                   'mhc.b',
                   'mhc.class',
                   'reference.id',
                   'species',
                   'v.segm',
                   'vdjdb.score']


def generate_slim_db(default_db: pd.DataFrame):
    slim_db = default_db[slim_db_columns]
    slim_db['j.start'] = default_db['cdr3fix'].apply(lambda x: x['jStart'] if not x == '' else x)
    slim_db['v.end'] = default_db['cdr3fix'].apply(lambda x: x['vEnd'] if not x == '' else x)
    slim_db.dropna(subset=slim_db_columns, inplace=True)
    slim_db.to_csv('../database/vdjdb_slim.txt', sep='\t')
