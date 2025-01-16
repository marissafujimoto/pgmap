from Levenshtein import distance, hamming
from Bio import Align

_blast_aligner = Align.PairwiseAligner()
_blast_aligner.match_score = 1.0
_blast_aligner.mismatch_score = -2.0
_blast_aligner.gap_score = -2.5


def hamming_score(candidate: str, reference: str) -> int:
    """
    Calculate the hamming score between two sequences. The hamming score is defined as the total length minus the
    number of mismatched base pairs.

    Args:
        candidate (str): The candidate sequence.
        reference (str): The reference sequence.

    Returns:
        int: The hamming score between the candidate and reference.
    """
    return len(candidate) - hamming(candidate, reference)


def edit_distance_score(candidate: str, reference: str) -> int:
    """
    Calculate the edit distance score between two sequences. The hamming score is defined as the total length minus
    the minimum number of insertions, deletions, and subsitutions to change the candidate sequence into the reference.

    Args:
        candidate (str): The candidate sequence.
        reference (str): The reference sequence.

    Returns:
        int: The edit distance score between the candidate and reference.
    """
    return len(reference) - distance(candidate, reference)


def blast_aligner_score(candidate: str, reference: str) -> int:
    """
    Calculate the blast alignment score between two sequences.

    Args:
        candidate (str): The candidate sequence.
        reference (str): The reference sequence.

    Returns:
        int: The edit distance score between the candidate and reference.
    """
    # TODO remove this? It might be useful as a library function, but is usused in the counts algorithm
    return _blast_aligner.score(candidate, reference)
