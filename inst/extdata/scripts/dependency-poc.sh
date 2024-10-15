#!/bin/bash

# TODO move source to /usr/src and build binaries to /usr/local/bin?

/home/seqtk/seqtk-1.4/seqtk trimfq /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz > /tmp/trimmed.fastq && echo "seqtk installed"

/home/config/idemp/idemp -b /home/config/barcode_ref_file.sample.txt -I1 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -R1 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -R2 /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -o /tmp/demux/ && echo "idemp installed"

# TODO pgPEN_index not compatible with bowtie2. Can we generate them?
# /home/biodocker/bin/bowtie2 -x /home/config/pgPEN_index/pgPEN_index/pgPEN_gRNA1.1.ebwt -q /home/data/tutorial_fastqs/PP_pgPEN_HeLa_S1_R1_001.fastq.gz -S /tmp/gRNA1.sam > /tmp/alignment.log && echo "bowtie2 installed"

/home/biodocker/bin/bowtie2 --version && echo "bowtie2 Installed"

/usr/src/samtools-1.9/samtools && echo "samtools Installed"
