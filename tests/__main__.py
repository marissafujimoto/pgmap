import unittest

from pgmap.io import fastq_reader, barcode_reader
from pgmap.trimming import read_trimmer


# TODO separate into different classes for modules
class TestIO(unittest.TestCase):

    def test_read_fastq(self):
        count = 0

        # TODO move these files to a cleaner test data folder
        for sequence in fastq_reader.read_fastq("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz"):
            count += 1

        self.assertEqual(count, 10000)

    def test_read_fastq_gz(self):
        count = 0

        for sequence in fastq_reader.read_fastq("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz"):
            count += 1

        self.assertEqual(count, 10000)

    def test_three_read_trim(self):
        count = 0

        for paired_read in read_trimmer.three_read_trim("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq.gz",
                                                        "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz",
                                                        "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I2_001_Sampled10k.fastq.gz"):
            count += 1

            # TODO test content of paired_read when code to read in library and barcodes is ready. For now manually inspect
            if count < 10:
                print(paired_read)

        self.assertEqual(count, 10000)

    def test_two_read_trim(self):
        count = 0

        for paired_read in read_trimmer.two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                                      "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"):
            count += 1

            # TODO test content of paired_read when code to read in library and barcodes is ready. For now manually inspect
            if count < 10:
                print(paired_read)

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


if __name__ == "__main__":
    unittest.main()
