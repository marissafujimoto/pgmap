#!/bin/bash

seqtk trimfq /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz > /tmp/trimmed.fastq && echo "seqtk installed"

idemp -b /home/config/barcode_ref_file.sample.txt -I1 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -R1 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -R2 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -o /tmp/demux/ && echo "idemp installed"

# Generating indices
# TODO: are these the right files for making indices for the samples here?
bowtie2-build /home/data/ref/pgPEN_R1.fa,/home/data/ref/pgPEN_R2.fa pg_PEN_index

bowtie2 -x pg_PEN_index \
  -1 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz \
  -2 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R2_001.fastq.gz \
  -S /tmp/gRNA1.sam > /tmp/alignment.log && echo "bowtie2 installed"

samtools && echo "samtools Installed"
