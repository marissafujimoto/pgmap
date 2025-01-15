from typing import Iterable


def construct_grna_error_alignment_cache(gRNAs: Iterable[str], gRNA_error_tolerance: int) -> dict[str, tuple[str, int]]:
    """
    Construct an alignment cache object containing all gRNAs with the error tolerance amount of substitutions. The 
    number of gRNAs within the error tolerance grows exponentially. As such this function should only be used for
    error tolerances of 0 or 1.

    Args:
        gRNAs (Iterable[str]): An iterable producing all the gRNAs to construct the alignment cache from.
        gRNA_error_tolerance (int): The error tolerance used to create gRNA alignment candidates.

    Returns:
        alignment_cache (dict[str, tuple[str, int]]): A mapping from each valid alignment string to a tuple of the
        reference alignment sequence and the hamming distance from the reference. Guarantees a minimum alignment,
        though there could be multiple best alignments depending on the gRNAs.

    Raises:
        ValueError: if gRNA_error_tolerance is not 0 or 1.
    """

    if gRNA_error_tolerance > 1 or gRNA_error_tolerance < 0:
        raise ValueError(f"gRNA error tolerance must be 0 or 1 but was {
                         gRNA_error_tolerance}")

    alignment_cache = {}

    if gRNA_error_tolerance:
        for gRNA in gRNAs:
            alignment_cache.update(
                {mutation: (gRNA, 1) for mutation in _get_mutations(gRNA)})

    alignment_cache.update({gRNA: (gRNA, 0) for gRNA in gRNAs}
                           )  # prefer perfect alignments

    return alignment_cache


def _get_mutations(gRNA: str) -> Iterable[str]:
    for i in range(len(gRNA)):
        for substitution in "ATCG":
            if gRNA[i] == substitution:
                continue

            yield gRNA[:i] + substitution + gRNA[i + 1:]
