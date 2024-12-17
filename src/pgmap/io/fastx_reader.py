from Bio import SeqIO
from mimetypes import guess_type

import gzip
from typing import IO, Iterable


def read_fastq(fastq_path: str) -> Iterable[str]:
    """
    Read sequences from a fastq (ignoring quality). Streams from a file and uses O(1) memory.

    Args:
        fastq_path (str): The path to a fastq file. Can optionally be gzipped.

    Yields:
        str: The next sequence in the fastq file.
    """
    # TODO check file validity?
    yield from _read_fastx(fastq_path, "fastq")


def read_fasta(fasta_path: str) -> Iterable[str]:
    """
    Read sequences from a fasta. Streams from a file and uses O(1) memory.

    Args:
        fasta_path (str): The path to a fasta file. Can optionally be gzipped.

    Yields:
        str: The next sequence in the fasta file.
    """
    # TODO check file validity?
    yield from _read_fastx(fasta_path, "fasta")


def _read_fastx(fastx_path: str, fastx_format_name: str) -> Iterable[str]:
    encoding = guess_type(fastx_path)[1]

    if encoding == "gzip":
        with gzip.open(fastx_path, "rt") as fastx_file:
            yield from _generate_sequences(fastx_file, fastx_format_name)
    else:
        with open(fastx_path) as fastx_file:
            yield from _generate_sequences(fastx_file, fastx_format_name)


def _generate_sequences(fastq_file: IO, fastx_format_name: str) -> Iterable[str]:
    for read in SeqIO.parse(fastq_file, fastx_format_name):
        yield str(read.seq)
