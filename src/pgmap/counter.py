from typing import Counter, Iterable

from pgmap.model.paired_read import PairedRead


def get_counts(candidates: Iterable[PairedRead], gRNA1_references: set[str], gRNA2_references: set[str], barcode_references: set[str]) -> Counter[tuple[str, str, str]]:
    # TODO implement counting and alignment algorithm
    raise NotImplementedError("Counts not implemented")
