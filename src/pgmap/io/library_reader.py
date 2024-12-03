from pgmap.io import fastx_reader

def read_paired_guide_library(R1_path: str, R2_path: str) -> tuple[set[str], set[str]]:
    gRNA1s = set(fastx_reader.read_fasta(R1_path))
    gRNA2s = set(fastx_reader.read_fasta(R2_path))

    return gRNA1s, gRNA2s