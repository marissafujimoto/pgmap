import unittest

from pgmap.io import fastq_reader


class TestIO(unittest.TestCase):

    def test_read_fastq(self):
        count = 0

        # TODO move these files to a cleaner test data folder
        for sequence in fastq_reader.read_fastq("inst/extdata/data/fastqs/fastq_trimmed/pgMAP_tutorial_gRNA1_trimmed.fastq"):
            count += 1

        self.assertEqual(count, 10e5)

    def test_read_fastq_gz(self):
        count = 0

        for sequence in fastq_reader.read_fastq("inst/extdata/data/fastqs/PP_pgPEN_HeLa_S1_R3_001.fastq.gz"):
            count += 1

        self.assertEqual(count, 10e5)


if __name__ == "__main__":
    unittest.main()
