from OneSideFixerResult import OneSideFixerResult
import FixType
from KmerScanner import KmerScanner
from FixerResults import FixerResult
from Utils import translate_linear, simplify_segment_name

import pandas as pd

from collections import defaultdict
from typing import Tuple, Optional, Any


class Cdr3Fixer:
    def __init__(self, segments_file_name: str, segments_seq_file_name: str,
                 max_replace_size: int = 1, min_hit_size: int = 2):
        self.segments_by_id_by_species = defaultdict(dict)
        self.segments_by_sequence_part_by_species_gene = defaultdict(dict)
        self.max_replace_size = max_replace_size
        self.min_hit_size = min_hit_size

        self._load_segments_data(segments_file_name)
        self._load_segments_sequence_data(segments_seq_file_name)
        self.nomenclature_conversions = defaultdict(lambda: None,
                                                    pd.read_csv("../patches/nomenclature.conversions",
                                                                sep='\t',
                                                                index_col=0,
                                                                header=None,
                                                                skiprows=1
                                                                )[1].to_dict()) #rewrite it

    def _load_segments_data(self, segments_file_name: str) -> None:

        segments_file = pd.read_csv(segments_file_name, sep='\t')
        for columns in ['#species', 'segment']:
            segments_file[columns] = segments_file[columns].apply(lambda x: x.lower())

        for _, segment in segments_file.iterrows():

            if segment.segment.lower().startswith("v") or segment.segment.lower().startswith("j"):
                is_j_segment = segment.segment.lower().startswith("j")
                segment.sequence = segment.sequence[:segment.reference_point + 4] if is_j_segment \
                    else segment.sequence[segment.reference_point - 3:]

                self.segments_by_id_by_species[segment['#species']][segment.id] = translate_linear(
                    segment.sequence, is_j_segment)

    def _load_segments_sequence_data(self, segments_seq_file_name: str) -> None:

        segments_seq_file = pd.read_csv(segments_seq_file_name, sep='\t')
        for _, segment in segments_seq_file.iterrows():
            species_chain = segment.species + (".alpha" if segment.gene == "TRA" else ".beta")
            self.segments_by_sequence_part_by_species_gene[species_chain][segment.cdr3] = segment.segm

    def get_closest_id(self, species: str, segment_id: str) -> Optional[str]:
        segments_by_id = self.segments_by_id_by_species[species.lower()]

        if species.lower() == "homosapiens":
            conversion = self.nomenclature_conversions[segment_id]
            if conversion:
                segment_id = conversion

        if not segments_by_id:
            return None

        for id_variant in [segment_id, *simplify_segment_name(segment_id)]:
            for possible_variant in [id_variant, f"{id_variant}*01", *[f"{id_variant}-{i}*01" for i in range(1, 101)]]:
                if possible_variant in segments_by_id.keys():
                    return possible_variant
        return None

    def get_segment_seq(self, species: str, segment_id: str) -> Optional[str]:
        segments_by_id = self.segments_by_id_by_species.get(species.lower())

        if not segments_by_id:
            return None

        return segments_by_id.get(segment_id)

    def fix(self, cdr3: str, segment_id: str, species: str, five_prime: bool) -> OneSideFixerResult:
        closest_id = self.get_closest_id(species, segment_id)
        segment_seq = self.get_segment_seq(species, closest_id)
        if not segment_seq:
            return OneSideFixerResult(cdr3 if five_prime else cdr3[::-1],
                                      closest_id,
                                      FixType.FailedBadSegment)
        if not five_prime:
            cdr3 = cdr3[::-1]
            segment_seq = segment_seq[::-1]

        scanner = KmerScanner(segment_seq, self.min_hit_size)

        hit = scanner.scan(cdr3)

        if hit:
            if hit.start_in_segment == 0:
                if hit.start_in_cdr3 == 0:
                    return OneSideFixerResult(cdr3, closest_id, FixType.NoFixNeeded, hit.match_size)
                else:
                    return OneSideFixerResult(cdr3[hit.start_in_cdr3:], closest_id, FixType.FixTrim, hit.match_size)
            else:
                if hit.start_in_cdr3 == 0:
                    return OneSideFixerResult(segment_seq[:hit.start_in_segment] + cdr3,
                                              closest_id,
                                              FixType.FixAdd,
                                              hit.match_size
                                              )
                elif hit.start_in_cdr3 <= self.max_replace_size:
                    return OneSideFixerResult(segment_seq[:hit.start_in_segment] + cdr3[hit.start_in_cdr3:],
                                              closest_id,
                                              FixType.FixReplace,
                                              hit.match_size
                                              )
                else:
                    return OneSideFixerResult(segment_seq,
                                              closest_id,
                                              FixType.FailedReplace,
                                              )
        else:
            return OneSideFixerResult(cdr3, closest_id, FixType.FailedNoAlignment)

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
                return ""
        else:
            for i in range(2, len(cdr3) - 3):
                res = segment_seq_map.get(cdr3[i:])
                if res:
                    return res
            else:
                return ""

    def fix_both(self, cdr3: str, v_id: str, j_id: str, species: str):
        v_results = [self.fix(cdr3, v, species, True) for v in v_id.split(",")]
        v_result = min(v_results, key=lambda x: x.FixType.rank)

        j_results = [self.fix(v_result.cdr3, j, species, False) for j in j_id.split(",")]
        j_result = min(j_results, key=lambda x: x.FixType.rank)

        new_cdr3 = j_result.cdr3[::-1]

        return FixerResult(
            new_cdr3,
            cdr3,
            new_cdr3 != cdr3,
            v_result.x,
            j_result.x if j_result.x < 0 else len(new_cdr3) - j_result.x,
            v_result.segmentId, v_result.FixType,
            j_result.segmentId, j_result.FixType
        )

