import pandas as pd

vdjdb_full = pd.read_csv('../database/vdjdb_full.txt', sep='\t',)
cluster_members = pd.read_csv('../database/cluster_members.txt', sep='\t')

vdjdb_full_clusters = vdjdb_full.copy()
vdjdb_full_clusters['cluster_member'] = 0

cluster_index_col = [
    'antigen.epitope',
    'species',
    'cdr3.alpha'
]

vdjdb_full_clusters = vdjdb_full_clusters.set_index(cluster_index_col)

cluster_members_tra = cluster_members[cluster_members['gene'] == 'TRA']

for _, cluster_member in cluster_members_tra.iterrows():
    vdjdb_full_clusters.loc[tuple(cluster_member[cluster_index_col[:-1] + ['cdr3aa']]), 'cluster_member'] = 1

vdjdb_full_clusters.reset_index(inplace=True)

cluster_index_col[-1] = 'cdr3.beta'

vdjdb_full_clusters = vdjdb_full_clusters.set_index(cluster_index_col)

cluster_members_trb = cluster_members[cluster_members['gene'] == 'TRB']

for _, cluster_member in cluster_members_trb.iterrows():
    vdjdb_full_clusters.loc[tuple(cluster_member[cluster_index_col[:-1] + ['cdr3aa']]), "cluster_member"] = 1

vdjdb_full_clusters = vdjdb_full_clusters.reset_index()[vdjdb_full.columns]
vdjdb_full_clusters.to_csv('../database/vdjdb_full_scored.txt', sep='\t')

slim_db = pd.read_csv('../database/vdjdb.slim.txt', sep='\t',)

cluster_index_col = [
    'antigen.epitope',
    'species',
    'gene'
]

slim_db_scored = slim_db.copy()
slim_db_scored = slim_db_scored.set_index(cluster_index_col)
slim_db_scored['cluster_member'] = 0

for _, cluster_member in cluster_members_tra.iterrows():
    slim_db_scored.loc[tuple(cluster_member[cluster_index_col]), 'cluster_member'] = 1

slim_db_scored = slim_db_scored.reset_index()[slim_db.columns]
slim_db_scored.to_csv('../database/vdjdb.slim.scored.txt', sep='\t',)