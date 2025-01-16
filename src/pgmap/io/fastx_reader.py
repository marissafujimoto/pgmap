from mimetypes import guess_type

import gzip
from typing import IO, Iterable
from contextlib import contextmanager


def read_fastq(fastq_path: str) -> Iterable[str]:
    """
    Read sequences from a fastq (ignoring quality). Streams from a file and uses O(1) memory.

    Args:
        fastq_path (str): The path to a fastq file. Can optionally be gzipped.

    Yields:
        str: The next sequence in the fastq file.
    """
    # TODO check file validity?
    with _gzip_agnostic_open(fastq_path) as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 1:
                yield line.strip()


def read_fasta(fasta_path: str) -> Iterable[str]:
    """
    Read sequences from a fasta. Streams from a file and uses O(1) memory.

    Args:
        fasta_path (str): The path to a fasta file. Can optionally be gzipped.

    Yields:
        str: The next sequence in the fasta file.
    """
    # TODO check file validity?
    with _gzip_agnostic_open(fasta_path) as fasta_file:
        for i, line in enumerate(fasta_file):
            if i % 2 == 1:
                yield line.strip()


@contextmanager
def _gzip_agnostic_open(path: str) -> Iterable[IO]:
    encoding = guess_type(path)[1]

    if encoding == "gzip":
        with gzip.open(path, "rt") as file:
            yield file
    else:
        with open(path) as file:
            yield file
