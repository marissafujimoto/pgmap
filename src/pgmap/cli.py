import argparse
import os
import sys

from pgmap.counter import counter
from pgmap.io import barcode_reader, library_reader, counts_writer
from pgmap.trimming import read_trimmer
from pgmap.model.trim_coordinate import TrimCoordinate
from pgmap.model.trim_strategy import TrimStrategy

TWO_READ_STRATEGY = "two-read"
THREE_READ_STRATEGY = "three-read"


def main():
    get_counts(_parse_args(sys.argv[1:]))


def get_counts(args: argparse.Namespace):
    barcodes = barcode_reader.read_barcodes(args.barcodes)
    gRNA1s, gRNA2s, gRNA_mappings, id_mapping = library_reader.read_paired_guide_library_annotation(
        args.library)

    candidate_reads = None

    candidate_reads = read_trimmer.trim(args.fastq, args.trim_strategy)

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
    parser.add_argument("--trim-strategy", required=True, type=_check_trim_strategy,
                        help="The trim strategy used to extract guides and barcodes. " +
                             "A custom trim strategy should be formatted as as space separate list of trim coordinates for gRNA1, gRNA2, and the barcode. " +
                             "Each trim coordinate should contain three zero indexed integers giving the file index relative to the order provided in --fastq, " +
                             "the inclusive start index of the trim, and the exclusive end index of the trim. " +
                             "For convenience the options \"two-read\" and \"three-read\" map to default values \"0:0:20 1:1:21 1:160:166\" and \"0:0:20 1:1:21 2:0:6\" respectively. " +
                             "The two read strategy should have fastqs R1 and I1. " +
                             "The three read strategy should have fastqs R1, I1, and I2.")
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


def _check_trim_strategy(serialized_trim_strategy: str) -> TrimStrategy:
    if serialized_trim_strategy == TWO_READ_STRATEGY:
        gRNA1_trim_coord = TrimCoordinate(file_index=0, start=0, end=20)
        gRNA2_trim_coord = TrimCoordinate(file_index=1, start=1, end=21)
        barcode_trim_coord = TrimCoordinate(file_index=1, start=160, end=166)

        return TrimStrategy(gRNA1=gRNA1_trim_coord, gRNA2=gRNA2_trim_coord, barcode=barcode_trim_coord)
    elif serialized_trim_strategy == THREE_READ_STRATEGY:
        gRNA1_trim_coord = TrimCoordinate(file_index=0, start=0, end=20)
        gRNA2_trim_coord = TrimCoordinate(file_index=1, start=1, end=21)
        barcode_trim_coord = TrimCoordinate(file_index=2, start=0, end=6)

        return TrimStrategy(gRNA1=gRNA1_trim_coord, gRNA2=gRNA2_trim_coord, barcode=barcode_trim_coord)
    else:
        serialized_trim_coordinates = serialized_trim_strategy.strip().split()

        if len(serialized_trim_coordinates) != 3:
            raise ValueError(
                f"Trim strategy {serialized_trim_strategy} must provide three coordinates but had {len(serialized_trim_coordinates)} coordinates instead")

        gRNA1_trim_coord, gRNA2_trim_coord, barcode_trim_coord = tuple(
            map(_check_trim_coordinate, serialized_trim_coordinates))

        return TrimStrategy(gRNA1=gRNA1_trim_coord, gRNA2=gRNA2_trim_coord, barcode=barcode_trim_coord)


def _check_trim_coordinate(serialized_trim_coordinate: str) -> TrimCoordinate:
    trim_indexes = serialized_trim_coordinate.split(":")

    if len(trim_indexes) != 3:
        raise ValueError(
            f"Trim coordinate {serialized_trim_coordinate} must provide three indexes but had {len(serialized_trim_coordinate)} indexes instead")

    file_index_str, start_str, end_str = trim_indexes

    if not file_index_str.isdigit() or not start_str.isdigit() or not end_str.isdigit():
        raise ValueError(
            f"Trim coordinate {serialized_trim_coordinate} contained a non integer trim index")

    file_index, start, end = tuple(map(int, trim_indexes))

    return TrimCoordinate(file_index=file_index, start=start, end=end)


if __name__ == "__main__":
    main()
