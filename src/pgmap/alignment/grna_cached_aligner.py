import itertools

from typing import Iterable


def construct_grna_error_alignment_cache(gRNAs: Iterable[str], gRNA_error_tolerance: int) -> dict[str, tuple[str, int]]:
    """
    Construct an alignment cache object containing all gRNAs with the error tolerance amount of substitutions. The 
    number of gRNAs within the error tolerance grows exponentially. As such this function should only be used for
    error tolerances of 0, 1, or 2.

    Args:
        gRNAs (Iterable[str]): An iterable producing all the gRNAs to construct the alignment cache from.
        gRNA_error_tolerance (int): The error tolerance used to create gRNA alignment candidates.

    Returns:
        alignment_cache (dict[str, tuple[str, int]]): A mapping from each valid alignment string to a tuple of the
        reference alignment sequence and the hamming distance from the reference. Guarantees a minimum alignment,
        though there could be multiple best alignments depending on the gRNAs.

    Raises:
        ValueError: if gRNA_error_tolerance is not 0, 1, or 2.
    """

    if gRNA_error_tolerance > 2 or gRNA_error_tolerance < 0:
        raise ValueError(
            "gRNA error tolerance must be 0, 1, or 2 but was " + str(gRNA_error_tolerance))

    alignment_cache = {}

    if gRNA_error_tolerance:
        for gRNA in gRNAs:
            # go from high subs to low subs to prefer better alignments
            for num_substitutions in reversed(range(1, gRNA_error_tolerance + 1)):
                alignment_cache.update(
                    {mutation: (gRNA, num_substitutions) for mutation in _get_mutations(gRNA, num_substitutions)})

    alignment_cache.update({gRNA: (gRNA, 0)
                           for gRNA in gRNAs})  # prefer perfect alignments

    return alignment_cache


def _get_mutations(gRNA: str, num_substitutions: int) -> Iterable[str]:
    for substitution_indices in _get_all_substitution_indices(gRNA, num_substitutions):
        yield from _generate_substitutions(gRNA, substitution_indices)


def _get_all_substitution_indices(gRNA: str, num_substitutions: int) -> Iterable[tuple[int]]:
    yield from itertools.combinations(range(len(gRNA)), num_substitutions)


def _generate_substitutions(gRNA: str, substitution_indices: Iterable[int]) -> Iterable[str]:
    for substitutions in itertools.product("ATCG", repeat=len(substitution_indices)):
        if any(gRNA[i] == substitution for i, substitution in zip(substitution_indices, substitutions)):
            continue

        bases = list(gRNA)

        for i, substitution in zip(substitution_indices, substitutions):
            bases[i] = substitution

        yield "".join(bases)
