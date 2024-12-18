from collections import defaultdict
import csv

from pgmap.io import fastx_reader


def read_paired_guide_library_fastas(R1_path: str, R2_path: str) -> tuple[set[str], set[str], dict[str, set[str]]]:
    """
    Reads a paired guide library from two fasta files. R1 and R2 are the first and second guide sequences
    respectively.

    Args:
        R1_path (str): The path to a fasta file for the first guide sequences.
        R2_path (str): The path to a fasta file for the second guide sequences.

    Returns:
        gRNA1s (set[str]): The set of all gRNA1 sequences.
        gRNA2s (set[str]): The set of all gRNA2 sequences.
        gRNA_mappings (dict[str, set[str]]): A mapping from each gRNA1 to a set of all of it's paired gRNA2s.
    """
    # TODO examine how other guide libraries are stored / shared. This might be niche
    # TODO deprecate?
    gRNA1s = set()
    gRNA2s = set()

    gRNA_mappings = defaultdict(set)

    for gRNA1, gRNA2 in zip(fastx_reader.read_fasta(R1_path), fastx_reader.read_fasta(R2_path)):
        gRNA1s.add(gRNA1)
        gRNA2s.add(gRNA2)

        gRNA_mappings[gRNA1].add(gRNA2)

    return gRNA1s, gRNA2s, gRNA_mappings


def read_paired_guide_library_annotation(annotation_path: str) -> tuple[set[str], set[str], dict[str, set[str]]]:
    """
    Reads a paired guide library from an annotation file. The file should be a tsv with a header of id, gRNA1, gRNA2.

    Args:
        annotation_path (str): The path to the library annotation file.

    Returns:
        gRNA1s (set[str]): The set of all gRNA1 sequences.
        gRNA2s (set[str]): The set of all gRNA2 sequences.
        gRNA_mappings (dict[str, set[str]]): A mapping from each gRNA1 to a set of all of it's paired gRNA2s.
    """
    # TODO docs

    gRNA1s = set()
    gRNA2s = set()

    gRNA_mappings = defaultdict(set)

    with open(annotation_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')

        for i, (id, gRNA1, gRNA2) in enumerate(tsv_reader):
            if i == 0:
                continue  # skip column headers

            gRNA1s.add(gRNA1)
            gRNA2s.add(gRNA2)

            gRNA_mappings[gRNA1].add(gRNA2)

    return gRNA1s, gRNA2s, gRNA_mappings
