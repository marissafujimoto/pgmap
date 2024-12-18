import argparse
import os
import sys

from pgmap.counter import counter

TWO_READ_STRATEGY = "two_read"
THREE_READ_STRATEGY = "three_read"


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
    parser.add_argument("--trim_strategy", required=False, choices=(TWO_READ_STRATEGY, THREE_READ_STRATEGY),
                        help="The trim strategy used to extract guides and barcodes. The two read strategy should have fastqs R1 and I1. The three read strategy should have fastqs R1, I1, and I2")  # TODO extract consts
    parser.add_argument("--gRNA2_error", required=False, default=2, type=_check_nonnegative,
                        help="The number of substituted base pairs to allow in gRNA2.")
    parser.add_argument("--barcode_error", required=False, default=2, type=_check_nonnegative,
                        help="The number of insertions, deletions, and subsititions of base pairs to allow in the barcodes.")
    return parser.parse_args(args)


def _check_nonnegative(value: str) -> int:
    int_value = int(value)

    if int_value < 0:
        raise ValueError(f"Count must be nonnegative but was {value}")

    return int_value


def _check_file_exists(path: str) -> str:
    if os.path.exists(path) and os.access(path, os.R_OK):
        return path
    else:
        raise argparse.ArgumentTypeError(
            f"File path {path} does not exist or is not readable")


if __name__ == "__main__":
    args = _parse_args(sys.argv[1:])
