import os
import argparse
import pandas as pd
import warnings
from termcolor import cprint


from ChunkQC import ChunkQC, ALL_COLS
from Cdr3Fixer import Cdr3Fixer
from DefaultDBGenerator import generate_default_db
from SlimDBGenerator import generate_slim_db

# reading and aggregating antigen nomenclature patches
antigen_df = pd.read_csv("../patches/antigen_epitope_species_gene.dict", sep="\t", index_col=0)
aggregated_species = antigen_df["antigen.species"].to_dict()
aggregated_gene = antigen_df["antigen.gene"].to_dict()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Arguments for database building")
    parser.add_argument("--chunks-to-build", help="chunks to include in database", nargs="+", type=str)
    parser.add_argument("--no2fix", action="store_true", help="Fix did not occurred if enabled")
    args = parser.parse_args()

    cprint("Reading, concatenating, and QC of chunks", "magenta")

    chunk_files = set(os.listdir("../chunks")).intersection(args.chunks_to_build) if \
        args.chunks_to_build else os.listdir("../chunks")

    chunk_files = [chunk_file for chunk_file in chunk_files if chunk_file[0] != "." and chunk_file.endswith(".txt")]

    cprint(f"Total number of chunks: {len(chunk_files)}", "magenta")
    chunk_df_list = []
    for chunk_file in chunk_files:
        chunk_df = pd.read_csv(f"../chunks/{chunk_file}", sep="\t", encoding_errors="ignore")
        chunk_qc = ChunkQC(chunk_df)
        chunk_error_messages = chunk_qc.process_chunk()

        if chunk_error_messages.keys():
            print(dict(chunk_error_messages))
            warn_message = f"There were errors processing {chunk_file}"
            warnings.warn(warn_message)

        chunk_df["antigen.species"] = chunk_df.T.apply(lambda x: aggregated_species.get(x["antigen.epitope"])
            if aggregated_species.get(x["antigen.epitope"]) else x["antigen.species"])
        chunk_df["antigen.gene"] = chunk_df.T.apply(lambda x: aggregated_gene.get(x["antigen.epitope"])
            if aggregated_gene.get(x["antigen.epitope"]) else x["antigen.gene"])


        chunk_df_list.append(chunk_df)

    os.makedirs("../database/", exist_ok=True)
    master_table = pd.concat(chunk_df_list)[ALL_COLS]

    cprint("Fixing CDR3 sequences (stage I)", "magenta")

    cdr3_fixer = Cdr3Fixer("../res/segments.txt", "../res/segments.aaparts.txt")

    for gene in ["alpha", "beta"]:
        master_table[f"v.{gene}"] = master_table.T.apply(
            lambda x: cdr3_fixer.guess_id(x[f"cdr3.{gene}"], x.species, gene, True
                                          ) if pd.isnull(x[f"v.{gene}"]) and not pd.isnull(x[f"cdr3.{gene}"])
            else x[f"v.{gene}"])

        master_table[f"j.{gene}"] = master_table.T.apply(
            lambda x: cdr3_fixer.guess_id(x[f"cdr3.{gene}"], x.species, gene, False
                                          ) if pd.isnull(x[f"j.{gene}"]) and not pd.isnull(x[f"cdr3.{gene}"])
            else x[f"j.{gene}"])

        fixer_results = master_table.T.apply(
            lambda x: cdr3_fixer.fix_both(x[f"cdr3.{gene}"],
                                          x[f"v.{gene}"],
                                          x[f"j.{gene}"],
                                          x.species,
                                          ) if not pd.isnull(x[f"cdr3.{gene}"]) else None)
        # remake fixer results
        master_table[f"cdr3.{gene}"] = fixer_results.apply(lambda x: x.cdr3 if x else None)
        master_table[f"v.{gene}"] = fixer_results.apply(lambda x: x.vId if x else None)
        master_table[f"j.{gene}"] = fixer_results.apply(lambda x: x.jId if x else None)
        master_table[f"cdr3fix.{gene}"] = fixer_results.apply(lambda x: x.results_to_dict() if x else None)

    master_table.set_index("cdr3.alpha").to_csv("../database/vdjdb_full.txt", sep="\t", quotechar='"')

    cprint("Generating and writing default database", "magenta")
    default_db = generate_default_db(master_table)

    cprint("Generating and writing slim database", "magenta")
    generate_slim_db(default_db)
    cprint("DB generation successfully finished!", "magenta")
