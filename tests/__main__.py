import unittest
import argparse
import csv
import os
import uuid
from collections import Counter

from pgmap.io import barcode_reader, fastx_reader, library_reader
from pgmap.trimming import read_trimmer
from pgmap.alignment import pairwise_aligner
from pgmap.counter import counter
from pgmap.model.paired_read import PairedRead
from pgmap import cli

THREE_READ_R1_PATH = "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq.gz"
THREE_READ_I1_PATH = "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz"
THREE_READ_I2_PATH = "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I2_001_Sampled10k.fastq.gz"
THREE_READ_BARCODES_PATH = "example-data/three-read-strategy/HeLa/screen_barcodes.txt"

TWO_READ_R1_PATH = "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz"
TWO_READ_I1_PATH = "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"
TWO_READ_BARCODES_PATH = "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/barcode_ref_file_revcomp.txt"

PGPEN_ANNOTATION_PATH = "example-data/pgPEN-library/paralog_pgRNA_annotations.txt"

INVALID_FILE_PATH = "file-does-not-exist"


class TestPgmap(unittest.TestCase):

    def test_read_fastq(self):
        count = 0

        for sequence in fastx_reader.read_fastq("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq"):
            count += 1

        self.assertEqual(count, 10000)

    def test_read_fastq_gz(self):
        count = 0

        for sequence in fastx_reader.read_fastq(THREE_READ_I1_PATH):
            count += 1

        self.assertEqual(count, 10000)

    def test_read_barcodes(self):
        barcodes = barcode_reader.read_barcodes(
            "example-data/three-read-strategy/HeLa/screen_barcodes.txt")

        self.assertEqual(len(barcodes), 6)
        self.assertEqual(barcodes["CTTGTA"], "sample1")
        self.assertEqual(barcodes["GGCTAC"], "sample2")
        self.assertEqual(barcodes["TAGCTT"], "sample3")
        self.assertEqual(barcodes["GATCAG"], "sample4")
        self.assertEqual(barcodes["ACTTGA"], "sample5")
        self.assertEqual(barcodes["GCCAAT"], "sample6")

    def test_read_library_fastas(self):
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library_fastas(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        self.assertEqual(len(gRNA1s), 5072)
        self.assertEqual(len(gRNA2s), 5095)

        # The sum of the length of all mappings should be equal to the number of reads in the library files
        self.assertEqual(sum(len(mapped_gRNA2s) for _, mapped_gRNA2s in gRNA_mappings.items()),
                         sum(1 for _ in fastx_reader.read_fasta("example-data/pgPEN-library/pgPEN_R1.fa")))

    def test_read_library_annotation(self):
        gRNA1s, gRNA2s, gRNA_mappings, id_mapping = library_reader.read_paired_guide_library_annotation(
            PGPEN_ANNOTATION_PATH)

        self.assertEqual(len(gRNA1s), 5072)
        self.assertEqual(len(gRNA2s), 5095)

        # The sum of the length of all mappings should be equal to the number of reads in the library files
        self.assertEqual(sum(len(mapped_gRNA2s) for _, mapped_gRNA2s in gRNA_mappings.items()),
                         sum(1 for _ in fastx_reader.read_fasta("example-data/pgPEN-library/pgPEN_R1.fa")))

        # The length of id_mapping should be equal to the line length of the file minus 1

        with open(PGPEN_ANNOTATION_PATH, 'r') as file:
            self.assertEqual(sum(1 for _ in csv.reader(
                file, delimiter='\t')) - 1, len(id_mapping))

        self.assertEqual(
            id_mapping[("ATTTCTATCCAAATCACTCA", "GAAAAAATTTGACTGCAGCA")], "AADAC_AADACL2_pg10")

    def test_three_read_trim(self):
        barcodes = barcode_reader.read_barcodes(THREE_READ_BARCODES_PATH)
        gRNA1s, gRNA2s, _ = library_reader.read_paired_guide_library_fastas(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.three_read_trim(THREE_READ_R1_PATH, THREE_READ_I1_PATH, THREE_READ_I2_PATH):
            count += 1

            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA2s and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        self.assertEqual(count, 10000)
        self.assertGreater(perfect_alignments / count, .7)

    def test_two_read_trim(self):
        barcodes = barcode_reader.read_barcodes(THREE_READ_BARCODES_PATH)
        gRNA1s, gRNA2s, _ = library_reader.read_paired_guide_library_fastas(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim(TWO_READ_R1_PATH, TWO_READ_I1_PATH):
            count += 1

            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA2s and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        self.assertEqual(count, 10000)
        self.assertGreater(perfect_alignments / count, .2)

    def test_hamming_score(self):
        self.assertEqual(pairwise_aligner.hamming_score("ABC", "AXX"), 1)
        self.assertEqual(pairwise_aligner.hamming_score("ABC", "XXX"), 0)
        self.assertEqual(pairwise_aligner.hamming_score("ABC", "ABC"), 3)

    def test_edit_distance_score(self):
        self.assertEqual(
            pairwise_aligner.edit_distance_score("ABCDEF", "XABCDE"), 4)
        self.assertEqual(
            pairwise_aligner.edit_distance_score("ABC", "XXX"), 0)
        self.assertEqual(
            pairwise_aligner.edit_distance_score("ABC", "ABC"), 3)

    def test_blast_aligner_score(self):
        self.assertEqual(
            pairwise_aligner.blast_aligner_score("ABC", "XAB"), -3)
        self.assertEqual(
            pairwise_aligner.blast_aligner_score("ABC", "XXX"), -6)
        self.assertEqual(
            pairwise_aligner.blast_aligner_score("ABC", "ABC"), 3)

    def test_counter_no_error_tolerance(self):
        barcodes = barcode_reader.read_barcodes(TWO_READ_BARCODES_PATH)
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library_fastas(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        candidate_reads = read_trimmer.two_read_trim(
            TWO_READ_R1_PATH, TWO_READ_I1_PATH)

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes, gRNA2_error_tolerance=0, barcode_error_tolerance=0)

        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim(TWO_READ_R1_PATH, TWO_READ_I1_PATH):
            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA_mappings[paired_read.gRNA1_candidate] and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        # Zero error tolerance should equal perfect alignment
        self.assertEqual(sum(paired_guide_counts.values()), perfect_alignments)

    def test_counter_default_error_tolerance(self):
        barcodes = barcode_reader.read_barcodes(TWO_READ_BARCODES_PATH)
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library_fastas(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        candidate_reads = read_trimmer.two_read_trim(
            TWO_READ_R1_PATH, TWO_READ_I1_PATH)

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes)

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim(TWO_READ_R1_PATH, TWO_READ_I1_PATH):
            count += 1

            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA_mappings[paired_read.gRNA1_candidate] and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        # Default error tolerance should be greater than perfect alignment but less than all counts
        self.assertGreater(
            sum(paired_guide_counts.values()), perfect_alignments)
        self.assertLess(
            sum(paired_guide_counts.values()), count)

    def test_counter_hardcoded_test_data(self):
        barcodes = {"COOL", "WOOD", "FOOD"}
        gRNA_mappings = {"LET": {"WOW", "LEG", "EAT"}}
        candidate_reads = [PairedRead("LET", "ROT", "FOOD"),
                           PairedRead("LET", "EXT", "FXOD"),
                           PairedRead("RUN", "LEG", "WOOD")]

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes)

        self.assertEqual(paired_guide_counts[("LET", "EAT", "FOOD")], 1)
        self.assertEqual(paired_guide_counts[("LET", "WOW", "FOOD")], 1)
        self.assertEqual(sum(paired_guide_counts.values()), 2)

    # TODO separate these into own test module?
    def test_arg_parse_happy_case(self):
        args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                "--library", PGPEN_ANNOTATION_PATH,
                                "--barcodes", TWO_READ_BARCODES_PATH,
                                "--trim_strategy", "two_read"])
        self.assertEqual(args.fastq[0], TWO_READ_R1_PATH)
        self.assertEqual(args.fastq[1], TWO_READ_I1_PATH)
        self.assertEqual(args.library, PGPEN_ANNOTATION_PATH)
        self.assertEqual(args.barcodes, TWO_READ_BARCODES_PATH)
        self.assertEqual(args.trim_strategy, "two_read")
        self.assertEqual(args.gRNA2_error, 2)
        self.assertEqual(args.barcode_error, 2)

    def test_arg_parse_invalid_fastq(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", INVALID_FILE_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two_read"])

    def test_arg_parse_invalid_library(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", INVALID_FILE_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two_read"])

    def test_arg_parse_invalid_barcodes(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", INVALID_FILE_PATH,
                                    "--trim_strategy", "two_read"])

    def test_arg_parse_invalid_trim_strategy(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "burger"])

    def test_arg_parse_negative_gRNA2_error(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two-read",
                                    "--gRNA2_error", "-1"])

    def test_arg_parse_negative_barcode_error(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two-read",
                                    "--barcode_error", "-1"])

    def test_arg_parse_invalid_type_gRNA2_error(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two-read",
                                    "--gRNA2_error", "one"])

    def test_arg_parse_invalid_type_barcode_error(self):
        with self.assertRaises(argparse.ArgumentError):
            args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                    "--library", PGPEN_ANNOTATION_PATH,
                                    "--barcodes", TWO_READ_BARCODES_PATH,
                                    "--trim_strategy", "two-read",
                                    "--barcode_error", "one"])

    def setUp(self):
        self.test_output_path = f"test_file_{uuid.uuid4().hex}.tsv"

    def tearDown(self):
        if os.path.exists(self.test_output_path):
            os.remove(self.test_output_path)

    def test_main_happy_case(self):
        args = cli._parse_args(["--fastq", TWO_READ_R1_PATH, TWO_READ_I1_PATH,
                                "--library", PGPEN_ANNOTATION_PATH,
                                "--barcodes", TWO_READ_BARCODES_PATH,
                                "--output", self.test_output_path,
                                "--trim_strategy", "two_read"])
        cli.get_counts(args)

        self.assertTrue(os.path.exists(self.test_output_path))

        barcodes = barcode_reader.read_barcodes(TWO_READ_BARCODES_PATH)

        gRNA1s, gRNA2s, gRNA_mappings, id_mapping = library_reader.read_paired_guide_library_annotation(
            PGPEN_ANNOTATION_PATH)

        candidate_reads = read_trimmer.two_read_trim(
            TWO_READ_R1_PATH, TWO_READ_I1_PATH)

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes)

        file_counts = self.load_counts_from_output_file(
            self.test_output_path, barcodes)

        self.assertEqual(file_counts, paired_guide_counts)

    def load_counts_from_output_file(self, output_path: str, barcodes: dict[str, str]) -> Counter[tuple[str, str, str]]:
        with open(output_path, 'r') as file:
            tsv_reader = csv.reader(file, delimiter='\t')

            sample_ids = None
            sample_id_to_barcode = {v: k for k, v in barcodes.items()}
            counts = Counter()

            for i, row in enumerate(tsv_reader):
                if i == 0:
                    sample_ids = row[3:]
                    continue

                gRNA1 = row[1]
                gRNA2 = row[2]

                for j in range(3, len(row)):
                    barcode = sample_id_to_barcode[sample_ids[j - 3]]

                    count = int(row[j])

                    counts[(gRNA1, gRNA2, barcode)] += count

            return counts


if __name__ == "__main__":
    unittest.main()
