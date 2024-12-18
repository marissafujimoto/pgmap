from collections import Counter
from typing import Counter, Iterable

from pgmap.model.paired_read import PairedRead
from pgmap.alignment import pairwise_aligner


def get_counts(paired_reads: Iterable[PairedRead],
               gRNA_mappings: dict[str, set[str]],
               barcodes: set[str],
               gRNA2_error_tolerance: int = 2,
               barcode_error_tolerance: int = 2) -> Counter[tuple[str, str, str]]:
    """
    Count paired guides for each sample barcode with tolerance for errors. gRNA1 matchs only through
    perfect alignment. gRNA2 aligns if there is a match of a known (gRNA1, gRNA2) pairing having hamming distance
    within the gRNA2 error tolerance. Finally the barcode aligns if there is a match aligned by edit distance
    within a separate barcode error tolerance.

    Args:
        paired_reads (Iterable[PairedRead]): An iterable producing the candidate reads to be counted. Can be
        generator to minimize memory usage.
        gRNA_mappings (dict[str, set[str]]): The known mappings of each reference library gRNA1 to the set of gRNA2s
        the gRNA1 is paired with.
        barcodes (set[str]): The sample barcode sequences.
        gRNA2_error_tolerance (int): The error tolerance for the hamming distance a gRNA2 candidate can be to the
        reference gRNA2. Defaults to 2.
        barcode_error_tolerance (int): The error tolerance for the edit distance a barcode candidate can be to the
        reference barcode. Defaults to 2.

    Returns:
        paired_guide_counts (Counter[tuple[str, str, str]]): The counts of each (gRNA1, gRNA2, barcode) detected
        within the paired reads.
    """
    # TODO should this keep track of metrics for how many paired reads get discarded and at which step? maybe with a verbose logging option?
    # TODO should we key with sample id instead of barcode sequence?
    # TODO sanity test error tolerances on real data
    # TODO should alignment algorithm be user configurable?

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
