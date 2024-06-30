import pandas as pd

from collections import defaultdict
from typing import List, Tuple, Optional
import os


class Cdr3Fixer:
    def __init__(self, segments_file_name: str, segments_seq_file_name: str,
                 max_replace_size: int = 1, min_hit_size: int = 2):
        self.segments_by_id_by_species = defaultdict(dict)
        self.segments_by_sequence_part_by_species_gene = defaultdict(dict)
        self.max_replace_size = max_replace_size
        self.min_hit_size = min_hit_size

        self._load_segments_data(segments_file_name)
        self._load_segments_sequence_data(segments_seq_file_name)

    def _load_segments_data(self, segments_file_name: str) -> None:

        segments_file = pd.read_csv(segments_file_name, sep='\t')

        for _, segment in segments_file.iterrows():

            if segment.segment.lower().startswith("v") or segment.segment.lower().startswith("j"):
                is_j_segment = segment.segment.lower().startswith("v")
                segment.sequence = segment.sequence[:segment.reference_point + 4] if is_j_segment \
                    else segment.sequence[segment.reference_point - 3:]

                self.segments_by_id_by_species[segment['#species']][segment.id] = self.translate_linear(
                    segment.sequence, is_j_segment)

    def _load_segments_sequence_data(self, segments_seq_file_name: str) -> None:

        segments_seq_file_name = pd.read_csv(segments_seq_file_name, sep='\t')
        for _, segment in segments_seq_file_name.iterrows():
            species_chain = segment.species + (".alpha" if segment.gene == "TRA" else ".beta")
            self.segments_by_sequence_part_by_species_gene[species_chain][segment.cdr3] = segment.segm

    def get_closest_id(self, species: str, segment_id: str) -> Optional[str]:
        segments_by_id = self.segments_by_id_by_species.get(species.lower())

        if species.lower() == "homosapiens":
            conversion = self.nomenclature_conversions.get(segment_id)
            if conversion:
                segment_id = conversion

        if not segments_by_id:
            return None

        for id_variant in [segment_id, self.simplify_segment_name(segment_id)]:
            for it in [id_variant, f"{it}*01", *[f"{it}-{i}*01" for i in range(1, 101)]]:
                if it in segments_by_id:
                    return it

        return None

    def get_segment_seq(self, species: str, segment_id: str) -> Optional[str]:
        segments_by_id = self.segments_by_id_by_species.get(species.lower())

        if not segments_by_id:
            return None

        return segments_by_id.get(segment_id)

    def fix(self, cdr3: str, segment_id: str, species: str, five_prime: bool) -> Tuple[str, str, str]:
        closest_id = self.get_closest_id(species, segment_id)
        segment_seq = self.get_segment_seq(species, closest_id)

        if not segment_seq:
            return (
                cdr3 if five_prime else cdr3[::-1],
                closest_id,
                "FailedBadSegment"
            )

        if not five_prime:
            cdr3 = cdr3[::-1]
            segment_seq = segment_seq[::-1]

        hit = self.find_hit(cdr3, segment_seq)

        if hit:
            if hit[1] == 0:
                if hit[0] == 0:
                    return cdr3, closest_id, "NoFixNeeded"
                else:
                    return cdr3[hit[0]:], closest_id, "FixTrim"
            else:
                if hit[0] == 0:
                    return segment_seq[:hit[1]] + cdr3, closest_id, "FixAdd"
                elif hit[0] <= self.max_replace_size:
                    return segment_seq[:hit[1]] + cdr3[hit[0]:], closest_id, "FixReplace"
                else:
                    return cdr3, closest_id, "FailedReplace"
        else:
            return cdr3, closest_id, "FailedNoAlignment"

    def guess_id(self, cdr3: str, species: str, gene: str, five_prime: bool) -> str:
        species_gene = f"{species}.{gene}"
        segment_seq_map = self.segments_by_sequence_part_by_species_gene.get(species_gene)

        if not segment_seq_map:
            return ""

        if five_prime:
            for i in range(len(cdr3) - 4, 1, -1):
                res = segment_seq_map.get(cdr3[:i])
                if res:
                    return res
        else:
            for i in range(2, len(cdr3) - 3):
                res = segment_seq_map.get(cdr3[i:])
                if res:
                    return res

        return ""

    def fix_both(self, cdr3: str, v_id: str, j_id: str, species: str, gene: str) -> Tuple[str, str]:
        v_results = [self.fix(cdr3, v, species, True) for v in v_id.split(",")]
        v_result = min(v_results, key=lambda x: self.rank_fix_type(x[2]))

        j_results = [self.fix(v_result[0], j, species, False) for j in j_id.split(",")]
        j_result = min(j_results, key=lambda x: self.rank_fix_type(x[2]))

        new_cdr3 = j_result[0][::-1]

        return new_cdr3, v_result[0], j_result[0]

    @staticmethod
    def translate_linear(seq: str, is_j_segment: bool) -> str:
        # Simplified placeholder for translate_linear function
        return seq

    @staticmethod
    def simplify_segment_name(segment_id: str) -> str:
        # Placeholder for simplify_segment_name function
        return segment_id

    @staticmethod
    def find_hit(cdr3: str, segment_seq: str) -> Optional[Tuple[int, int]]:
        # Placeholder for find_hit function
        return None

    @staticmethod
    def rank_fix_type(fix_type: str) -> int:
        # Placeholder for rank_fix_type function
        return 0
