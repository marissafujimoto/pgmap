import unittest

from pgmap.io import barcode_reader, fastx_reader, library_reader
from pgmap.trimming import read_trimmer
from pgmap.alignment import pairwise_aligner
from pgmap.counter import counter
from pgmap.model.paired_read import PairedRead


class TestPgmap(unittest.TestCase):

    def test_read_fastq(self):
        count = 0

        for sequence in fastx_reader.read_fastq("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq"):
            count += 1

        self.assertEqual(count, 10000)

    def test_read_fastq_gz(self):
        count = 0

        for sequence in fastx_reader.read_fastq("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz"):
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

    def test_read_library(self):
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        self.assertEqual(len(gRNA1s), 5072)
        self.assertEqual(len(gRNA2s), 5095)

        # The sum of the lenght of all mappings should be equal to the number of reads in the library files
        self.assertEqual(sum(len(mapped_gRNA2s) for _, mapped_gRNA2s in gRNA_mappings.items()),
                         sum(1 for _ in fastx_reader.read_fasta("example-data/pgPEN-library/pgPEN_R1.fa")))

    def test_three_read_trim(self):
        barcodes = barcode_reader.read_barcodes(
            "example-data/three-read-strategy/HeLa/screen_barcodes.txt")
        gRNA1s, gRNA2s, _ = library_reader.read_paired_guide_library(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.three_read_trim("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq.gz",
                                                        "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz",
                                                        "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I2_001_Sampled10k.fastq.gz"):
            count += 1

            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA2s and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        self.assertEqual(count, 10000)
        self.assertGreater(perfect_alignments / count, .7)

    def test_two_read_trim(self):
        barcodes = barcode_reader.read_barcodes(
            "example-data/three-read-strategy/HeLa/screen_barcodes.txt")
        gRNA1s, gRNA2s, _ = library_reader.read_paired_guide_library(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                      "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"):
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
        barcodes = barcode_reader.read_barcodes(
            "example-data/three-read-strategy/HeLa/screen_barcodes.txt")
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        candidate_reads = read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                     "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz")

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes, gRNA2_error_tolerance=0, barcode_error_tolerance=0)

        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                      "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"):
            if paired_read.gRNA1_candidate in gRNA1s and paired_read.gRNA2_candidate in gRNA_mappings[paired_read.gRNA1_candidate] and paired_read.barcode_candidate in barcodes:
                perfect_alignments += 1

        # Zero error tolerance should equal perfect alignment
        self.assertEqual(sum(paired_guide_counts.values()), perfect_alignments)

    def test_counter_default_error_tolerance(self):
        barcodes = barcode_reader.read_barcodes(
            "example-data/three-read-strategy/HeLa/screen_barcodes.txt")
        gRNA1s, gRNA2s, gRNA_mappings = library_reader.read_paired_guide_library(
            "example-data/pgPEN-library/pgPEN_R1.fa", "example-data/pgPEN-library/pgPEN_R2.fa")

        candidate_reads = read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                     "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz")

        paired_guide_counts = counter.get_counts(
            candidate_reads, gRNA_mappings, barcodes)

        count = 0
        perfect_alignments = 0

        for paired_read in read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                      "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"):
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


if __name__ == "__main__":
    unittest.main()
