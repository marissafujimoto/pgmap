from collections import Counter
from typing import Counter, Iterable

from pgmap.model.paired_read import PairedRead
from pgmap.alignment import pairwise_aligner


def get_counts(paired_reads: Iterable[PairedRead],
               gRNA_mappings: dict[str, set[str]],
               barcodes: set[str],
               gRNA2_error_tolerance: int = 2,
               barcode_error_tolerance: int = 2) -> Counter[tuple[str, str, str]]:
    # TODO should this keep track of metrics for how many paired reads get discarded and at which step? maybe with a verbose logging option?

    paired_guide_counts = Counter()

    for paired_read in paired_reads:
        gRNA1 = paired_read.gRNA1_candidate

        if gRNA1 not in gRNA_mappings:
            continue

        gRNA2_score, gRNA2 = max((pairwise_aligner.hamming_score(paired_read.gRNA2_candidate, reference), reference)
                                 for reference in gRNA_mappings[paired_read.gRNA1_candidate])

        if (len(gRNA2) - gRNA2_score) > gRNA2_error_tolerance:
            continue

        barcode_score, barcode = max((pairwise_aligner.edit_distance_score(
            paired_read.barcode_candidate, reference), reference) for reference in barcodes)

        if (len(barcode) - barcode_score) > barcode_error_tolerance:
            continue

        # TODO data structure for this?
        paired_guide_counts[(gRNA1, gRNA2, barcode)] += 1

    return paired_guide_counts
