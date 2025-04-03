from Levenshtein import distance, hamming


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

