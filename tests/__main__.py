import unittest

from pgmap.io import fastq_reader
from pgmap.trimmer import two_read_trim, three_read_trim


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

        for paired_read in three_read_trim("example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_R1_001_Sampled10k.fastq.gz",
                                           "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I1_001_Sampled10k.fastq.gz",
                                           "example-data/three-read-strategy/HeLa/PP_pgRNA_HeLa_S1_I2_001_Sampled10k.fastq.gz"):
            count += 1

            # TODO test content of paired_read when code to read in library and barcodes is ready. For now manually inspect
            if count < 10:
                print(paired_read)

        self.assertEqual(count, 10000)

    def test_two_read_trim(self):
        count = 0

        for paired_read in two_read_trim("example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz",
                                         "example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz"):
            count += 1

            # TODO test content of paired_read when code to read in library and barcodes is ready. For now manually inspect
            if count < 10:
                print(paired_read)

        self.assertEqual(count, 10000)


if __name__ == "__main__":
    unittest.main()
