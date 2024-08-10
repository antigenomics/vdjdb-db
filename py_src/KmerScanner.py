from dataclasses import dataclass
from typing import Optional


@dataclass
class SearchResult:
    """
    Dataclass to store search results
    """
    start_in_segment: int
    start_in_cdr3: int
    match_size: int


class KmerScanner:
    """
    Scans other segment sequence to get overlap with cdr3 sequence
    """

    def __init__(self, seq, min_hit_size=2):
        """
        Initiate KmerScanner by creating set of kmers from segment sequence
        :param seq: amino acids sequence of the segment
        :param min_hit_size: minimal overlap between cdr3 sequence and segment sequence
        """
        self.min_hit_size = min_hit_size
        self.kmers = {}

        for i in range(min_hit_size, len(seq) + 1):
            for j in range(0, len(seq) - i + 1):
                kmer = seq[j:j + i]
                self.kmers[kmer] = j

    def scan(self, seq) -> Optional[SearchResult]:
        """
        scans other cdr3 sequences to get overlap with segment sequence
        :param seq: amino acids sequence of cdr3 to be fixed
        :return: search results with best overlap coordinates in cdr3
        """
        best_hit = None
        for i in range(self.min_hit_size, len(seq)):
            for j in range(len(seq) - i):
                kmer = seq[j:j + i]
                hit = self.kmers.get(kmer)
                if hit is not None:
                    current_hit = SearchResult(hit, j, len(kmer))
                else:
                    current_hit = SearchResult(-1, j, -1)

                if best_hit is None or (current_hit.match_size > best_hit.match_size):
                    best_hit = current_hit

        return best_hit if best_hit and best_hit.match_size > 0 else None
