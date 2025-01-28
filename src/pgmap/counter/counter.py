from collections import Counter
from typing import Counter, Iterable
import itertools

from pgmap.model.paired_read import PairedRead
from pgmap.alignment import pairwise_aligner, grna_cached_aligner


def get_counts(paired_reads: Iterable[PairedRead],
               gRNA_mappings: dict[str, set[str]],
               barcodes: set[str],
               gRNA1_error_tolerance: int = 1,
               gRNA2_error_tolerance: int = 1,
               barcode_error_tolerance: int = 1) -> Counter[tuple[str, str, str]]:
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

    gRNA1_cached_aligner = grna_cached_aligner.construct_grna_error_alignment_cache(
        gRNA_mappings.keys(), gRNA1_error_tolerance)

    gRNA2_cached_aligner = grna_cached_aligner.construct_grna_error_alignment_cache(
        set(itertools.chain.from_iterable(gRNA_mappings.values())), gRNA2_error_tolerance)

    for paired_read in paired_reads:
        paired_read.gRNA1_candidate

        if paired_read.gRNA1_candidate not in gRNA1_cached_aligner:
            continue

        gRNA1, _ = gRNA1_cached_aligner[paired_read.gRNA1_candidate]

        if paired_read.gRNA2_candidate not in gRNA2_cached_aligner:
            continue

        gRNA2, _ = gRNA2_cached_aligner[paired_read.gRNA2_candidate]

        if gRNA2 not in gRNA_mappings[gRNA1]:
            continue

        barcode_score, barcode = max((pairwise_aligner.edit_distance_score(
            paired_read.barcode_candidate, reference), reference) for reference in barcodes)

        if (len(barcode) - barcode_score) > barcode_error_tolerance:
            continue

        # TODO data structure for this?
        paired_guide_counts[(gRNA1, gRNA2, barcode)] += 1

    return paired_guide_counts
