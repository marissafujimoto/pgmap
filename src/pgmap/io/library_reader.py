from collections import defaultdict

from pgmap.io import fastx_reader


def read_paired_guide_library(R1_path: str, R2_path: str) -> tuple[set[str], set[str], dict[str, set[str]]]:
    gRNA1s = set()
    gRNA2s = set()

    gRNA_mappings = defaultdict(set)

    for gRNA1, gRNA2 in zip(fastx_reader.read_fasta(R1_path), fastx_reader.read_fasta(R2_path)):
        gRNA1s.add(gRNA1)
        gRNA2s.add(gRNA2)

        gRNA_mappings[gRNA1].add(gRNA2)

    return gRNA1s, gRNA2s, gRNA_mappings
