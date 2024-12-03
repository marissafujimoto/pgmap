from Bio import SeqIO
from mimetypes import guess_type

import gzip
from typing import IO, Iterable


def read_fastq(fastq_path: str) -> Iterable[str]:
    """
    Read sequences from a fastq (ignoring quality)

    Args:
        fastq_path (str): The path to a fastq file. Can optionally be gzipped

    Yields:
        str: The next sequence in the fastq file
    """
    # TODO check file validity?

    encoding = guess_type(fastq_path)[1]

    if encoding == "gzip":
        with gzip.open(fastq_path, "rt") as fastq_file:
            yield from _generate_sequences(fastq_file)
    else:
        with open(fastq_path) as fastq_file:
            yield from _generate_sequences(fastq_file)


def _generate_sequences(fastq_file: IO) -> Iterable[str]:
    for read in SeqIO.parse(fastq_file, "fastq"):
        yield str(read.seq)

# TODO add fasta reader for guide references