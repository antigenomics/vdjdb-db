import os
import argparse
import pandas as pd
import warnings
from dataclasses import asdict


from ChunkQC import ChunkQC, ALL_COLS
# from AlignBestSegments import *
from Cdr3Fixer import Cdr3Fixer
from GenerateDefaultDB import generate_default_db

antigen_df = pd.read_csv("../patches/antigen_epitope_species_gene.dict", sep='\t', index_col=0)
aggregated_species = antigen_df.groupby(level=0)['antigen.species'].agg(lambda x: x.iloc[0] if len(x) == 1 else list(x))
aggregated_gene = antigen_df.groupby(level=0)['antigen.gene'].agg(lambda x: x.iloc[0] if len(x) == 1 else list(x))

if __name__ == 'main':

    os.makedirs('../tmp/', exist_ok=True)

    parser = argparse.ArgumentParser(description='Arguments for database building')
    parser.add_argument('--chunks-to-build', help='chunks to include in database', nargs='+', type=str)
    parser.add_argument('--no2fix', action='store_true', help='Fix did not occurred if enabled')
    args = parser.parse_args()

    chunk_files = set(os.listdir('../chunks')).intersection(args.chunks_to_build) if \
        args.chunks_to_build else os.listdir('../chunks')

    chunk_files = [chunk_file for chunk_file in chunk_files if chunk_file[0] != '.' and chunk_file.endswith('.txt')]

    chunk_df_list = []

    for chunk_file in chunk_files:
        chunk_df = pd.read_csv(chunk_file, sep='\t', encoding_errors='ignore')
        chunk_qc = ChunkQC(chunk_df)
        chunk_error_messages = chunk_qc.process_chunk()

        if chunk_error_messages.keys():
            print(dict(chunk_error_messages))
            warn_message = f"There were errors processing {chunk_file}"
            warnings.warn(warn_message)

        chunk_df['antigen.species'] = chunk_df['antigen.epitope'].apply(
            lambda x: aggregated_species[x] if x in aggregated_species else None
        )
        chunk_df['antigen.gene'] = chunk_df['antigen.epitope'].apply(
            lambda x: aggregated_gene[x] if x in aggregated_gene else None
        )

        chunk_df_list.append(chunk_df)

    os.makedirs('../database/', exist_ok=True)
    master_table = pd.concat(chunk_df_list)[ALL_COLS]
    master_table.to_csv('../database/vdjdb_full.txt', sep='\t')

    cdr3_fixer = Cdr3Fixer("../res/segments.txt", "../res/segments.aaparts.txt")
    print("Fixing CDR3 sequences (stage I)")

    for gene in ['alpha', 'beta']:
        master_table[f'v.{gene}'] = master_table.T.apply(
            lambda x: cdr3_fixer.guess_id(x[f'cdr3.{gene}'], x.species, gene, True
                                          ) if pd.isnull(x[f'v.{gene}']) and not pd.isnull(x[f'cdr3.{gene}'])
            else x[f'v.{gene}'])

        master_table[f'j.{gene}'] = master_table.T.apply(
            lambda x: cdr3_fixer.guess_id(x[f'cdr3.{gene}'], x.species, gene, False
                                          ) if pd.isnull(x[f'j.{gene}']) and not pd.isnull(x[f'cdr3.{gene}'])
            else x[f'j.{gene}'])

        fixer_results = master_table.T.apply(
            lambda x: cdr3_fixer.fix(x[f'cdr3.{gene}'],
                                     x[f'v.{gene}'],
                                     x[f'j.{gene}'],
                                     x.species,
                                     ) if not pd.isnull(x[f'cdr3.{gene}']) else None)

        master_table[f'cdr3.{gene}'] = fixer_results.apply(lambda x: x.cdr3)
        master_table[f'cdr3fix.{gene}'] = fixer_results.apply(lambda x: asdict(x) if x else None)


    default_db = generate_default_db(master_table)
    print("Generating and writing slim database")
