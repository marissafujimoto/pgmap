from typing import Iterable


def construct_grna_error_alignment_cache(gRNAs: list[str], gRNA_error_tolerance: int) -> dict[str, tuple[str, int]]:
    # TODO docs

    # TODO throw error if tolerance > 1 for now
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
