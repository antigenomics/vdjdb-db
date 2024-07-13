from dataclasses import dataclass


@dataclass
class SearchResult:
    start_in_segment: int
    start_in_cdr3: int
    match_size: int

class KmerScanner:
    """
    Scans other segment seq to get overlap with cdr3 seq
    """

    def __init__(self, seq, min_hit_size=2):
        self.min_hit_size = min_hit_size
        self.kmers = {}

        for i in range(min_hit_size, len(seq) + 1):
            for j in range(0, len(seq) - i + 1):
                kmer = seq[j:j + i]
                self.kmers[kmer] = j

    def scan(self, seq) -> SearchResult:
        best_hit = max(
            (SearchResult(self.kmers.get(seq[j:j + i], -1), j, i) for i in range(self.min_hit_size, len(seq))
             for j in range(len(seq) - i + 1))
            , key=lambda x: x.match_size)

        return best_hit if best_hit.match_size > 0 else None

