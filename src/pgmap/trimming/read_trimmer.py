from typing import Iterable

from pgmap.io.fastx_reader import read_fastq
from pgmap.model.trim_coordinate import TrimCoordinate
from pgmap.model.trim_strategy import TrimStrategy
from pgmap.model.paired_read import PairedRead


def two_read_trim(R1_path: str, I1_path: str, ) -> Iterable[PairedRead]:
    """
    Trim raw fastqs into PairedReads using a hardcoded two read strategy. Streams from a file and uses O(1) memory.

    Attributes:
        fastq_paths (list[str]): A list of paths to raw fastq files to trim from.
        trim_strategy (TrimStrategy): A TrimStrategy object describing how to trim the raw reads.

    Yields:
        candidates (PairedRead): Trimmed paired reads.
    """
    # TODO include some kind of citation or link to the two read sequencing strategy?
    R1_trim_coord = TrimCoordinate(file_index=0, start=0, end=20)
    R2_trim_coord = TrimCoordinate(file_index=1, start=1, end=21)
    barcode_trim_coord = TrimCoordinate(file_index=1, start=160, end=166)

    trim_strategy = TrimStrategy(
        R1=R1_trim_coord, R2=R2_trim_coord, barcode=barcode_trim_coord)

    yield from trim([R1_path, I1_path], trim_strategy)


def three_read_trim(R1_path: str, R2_path: str, I1_path: str) -> Iterable[PairedRead]:
    """
    Trim raw fastqs into PairedReads using a hardcoded three read strategy. Streams from a file and uses O(1) memory.

    Attributes:
        fastq_paths (list[str]): A list of paths to raw fastq files to trim from.
        trim_strategy (TrimStrategy): A TrimStrategy object describing how to trim the raw reads.

    Yields:
        candidates (PairedRead): Trimmed paired reads.
    """
    # TODO include some kind of citation or link to the three read sequencing strategy?
    R1_trim_coord = TrimCoordinate(file_index=0, start=0, end=20)
    R2_trim_coord = TrimCoordinate(file_index=1, start=1, end=21)
    barcode_trim_coord = TrimCoordinate(file_index=2, start=0, end=6)

    trim_strategy = TrimStrategy(
        R1=R1_trim_coord, R2=R2_trim_coord, barcode=barcode_trim_coord)

    yield from trim([R1_path, R2_path, I1_path], trim_strategy)


def trim(fastq_paths: list[str], trim_strategy: TrimStrategy) -> Iterable[PairedRead]:
    """
    Trim raw fastqs into PairedReads. Uses O(1) memory. Streams from a file and uses O(1) memory.

    Attributes:
        fastq_paths (list[str]): A list of paths to raw fastq files to trim from.
        trim_strategy (TrimStrategy): A TrimStrategy object describing how to trim the raw reads.

    Yields:
        candidates (PairedRead): Trimmed paired reads.
    """
    for raw in zip(*list(map(read_fastq, fastq_paths))):
        candidate = PairedRead(gRNA1_candidate=raw[trim_strategy.R1.file_index][trim_strategy.R1.start:trim_strategy.R1.end],
                               gRNA2_candidate=raw[trim_strategy.R2.file_index][trim_strategy.R2.start:trim_strategy.R2.end],
                               barcode_candidate=raw[trim_strategy.barcode.file_index][trim_strategy.barcode.start:trim_strategy.barcode.end])

        yield candidate
