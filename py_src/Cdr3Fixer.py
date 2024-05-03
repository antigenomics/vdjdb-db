import pandas as pd


class Cdr3Fixer():
    def __init__(self, segments_file_name, segments_seq_file_name, max_replace_size=1, min_hit_size=2):
        self.segments_by_sequence_part_by_species_gene = {}
        self.segments_by_id_by_species = {}
        self.max_replace_size = max_replace_size
        self.min_hit_size = min_hit_size
        self.nomenclature_conversions = pd.read_csv("../patches/nomenclature.conversions",
                                                    sep='\t', index_col=0, header=False).set_index(0).to_dict()






