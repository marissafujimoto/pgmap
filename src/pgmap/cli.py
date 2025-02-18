import argparse
import os
import sys

from pgmap.counter import counter
from pgmap.io import barcode_reader, library_reader, counts_writer
from pgmap.trimming import read_trimmer

TWO_READ_STRATEGY = "two-read"
THREE_READ_STRATEGY = "three-read"


def main():
    get_counts(_parse_args(sys.argv[1:]))


def get_counts(args: argparse.Namespace):
    barcodes = barcode_reader.read_barcodes(args.barcodes)
    gRNA1s, gRNA2s, gRNA_mappings, id_mapping = library_reader.read_paired_guide_library_annotation(
        args.library)

    candidate_reads = None

    if args.trim_strategy == TWO_READ_STRATEGY:
        candidate_reads = read_trimmer.two_read_trim(*args.fastq)
    elif args.trim_strategy == THREE_READ_STRATEGY:
        candidate_reads = read_trimmer.two_read_trim(*args.fastq)

    paired_guide_counts = counter.get_counts(
        candidate_reads, gRNA_mappings, barcodes, gRNA1_error_tolerance=args.gRNA1_error, gRNA2_error_tolerance=args.gRNA2_error, barcode_error_tolerance=args.barcode_error)

    counts_writer.write_counts(
        args.output, paired_guide_counts, barcodes, id_mapping)


def _parse_args(args: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="pgmap", description="A tool to count paired guides from CRISPR double knockout screens.", exit_on_error=False)
    # TODO in general these file formats should be documented more
    parser.add_argument("-f", "--fastq", nargs='+', required=True, type=_check_file_exists,
                        help="Fastq files to count from, separated by space. Can optionally be gzipped.")
    parser.add_argument("-l", "--library", required=True, type=_check_file_exists,
                        help="File containing annotated pgRNA information including the pgRNA id and both guide sequences.")
    # TODO support no barcodes?
    parser.add_argument("-b", "--barcodes", required=True, type=_check_file_exists,
                        help="File containing sample barcodes including the barcode sequence and the sample id.")
    # TODO check can write to this path?
    parser.add_argument("-o", "--output", required=False,
                        help="Output file path to populate with the counts for each paired guide and sample. If not provided the counts will be output in STDOUT.")
    # TODO support arbitrary trim strategies
    parser.add_argument("--trim-strategy", required=True, choices=(TWO_READ_STRATEGY, THREE_READ_STRATEGY),
                        help="The trim strategy used to extract guides and barcodes. The two read strategy should have fastqs R1 and I1. The three read strategy should have fastqs R1, I1, and I2")  # TODO extract consts
    parser.add_argument("--gRNA1-error", required=False, default=1, type=_check_gRNA1_error,
                        help="The number of substituted base pairs to allow in gRNA1. Must be less than 3.")
    parser.add_argument("--gRNA2-error", required=False, default=1, type=_check_gRNA2_error,
                        help="The number of substituted base pairs to allow in gRNA2. Must be less than 3.")
    parser.add_argument("--barcode-error", required=False, default=1, type=_check_barcode_error,
                        help="The number of insertions, deletions, and subsititions of base pairs to allow in the barcodes.")
    return parser.parse_args(args)


def _check_gRNA1_error(value: str) -> int:
    int_value = int(value)

    if int_value < 0:
        raise ValueError(f"gRNA1-error must be nonnegative but was {value}")

    if int_value > 2:
        raise ValueError(f"gRNA1-error must be less than 3 but was {value}")

    return int_value


def _check_gRNA2_error(value: str) -> int:
    int_value = int(value)

    if int_value < 0:
        raise ValueError(f"gRNA2-error must be nonnegative but was {value}")

    if int_value > 2:
        raise ValueError(f"gRNA2-error must be less than 3 but was {value}")

    return int_value


def _check_barcode_error(value: str) -> int:
    int_value = int(value)

    if int_value < 0:
        raise ValueError(f"barcode-error must be nonnegative but was {value}")

    return int_value


def _check_file_exists(path: str) -> str:
    if os.path.exists(path) and os.access(path, os.R_OK):
        return path
    else:
        raise argparse.ArgumentTypeError(
            f"File path {path} does not exist or is not readable")


if __name__ == "__main__":
    main()
