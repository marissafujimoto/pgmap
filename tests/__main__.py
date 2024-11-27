import unittest

from pgmap.io import fastq_reader

class TestIO(unittest.TestCase):

    def test_hello_world(self):
        self.assertEqual(fastq_reader.hello_world(), "Hello World")

if __name__ == "__main__":
    unittest.main()
