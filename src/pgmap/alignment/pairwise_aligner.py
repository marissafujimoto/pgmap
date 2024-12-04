import editdistance
from Bio import Align

_blast_aligner = Align.PairwiseAligner()
_blast_aligner.match_score = 1.0
_blast_aligner.mismatch_score = -2.0
_blast_aligner.gap_score = -2.5


def hamming_score(candidate: str, reference: str) -> int:
    return sum(map(str.__eq__, candidate, reference))


def edit_distance_score(candidate: str, reference: str) -> int:
    return len(reference) - editdistance.eval(candidate, reference)


def blast_aligner_score(candidate: str, reference: str) -> int:
    return _blast_aligner.score(candidate, reference)
