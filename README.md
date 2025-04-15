# pgmap - (paired guide RNA MAPper) Pipeline

pgmap is a command line tool which processes raw sequencing data into tables of counts per genetic perturbation. pgmap is purpose built for genetic interaction analysis; rather than utilizing bioinformatics tools designed for general purpose alignment of DNA sequences, pgmap exploits the paired structure of genetic interaction to efficiently map and count genetic perturbations.

It is based off of the original code and research from the Berger Lab stored in this repository: [https://github.com/FredHutch/pgMAP_pipeline](https://github.com/FredHutch/pgMAP_pipeline)

## Prerequisites

pgmap requires python3 above or equal to version 3.10.

### Installation using pip

```
pip install pgmap
```

Check if pgmap is installed in path
```
pgmap --help
```

### Installation from Source Code

For development purposes you can also install pgmap from the source code.

First you can clone this repo
```
git clone https://github.com/FredHutch/pgmap
```

To install required dependencies
```
pip install -r requirements.txt
```

Now you can install the package
```
pip install .
```

Run tests to verify pgmap is running
```
python3 -m tests
```

Check if pgmap is installed in path:
```
pgmap --help
```

## Getting Started

You'll need to declare four files (we have example data to work with):

- Two sequencing reads files:
  - I1 fastq file
  - R1 fastq file
- Barcodes path which says which barcodes go to which cells
- pgPEN annotations path (available in `example-data`) OR your own library annotations (see [Library Annotations](#library-annotations) for the format).
- A trim strategy. If using pgPEN and a two read strategy, use `two-read`. Otherwise see [Custom Trim Strategies](#custom-trim-strategies).

## pgmap CLI

| Argument | Description |
| - | - |
| -f FASTQ [FASTQ ...], --fastq FASTQ [FASTQ ...] | Fastq files to count from, separated by space. Can optionally be gzipped. |
| -l LIBRARY, --library LIBRARY | File containing annotated pgRNA information including the pgRNA id and both guide sequences. |
| -b BARCODES, --barcodes BARCODES | File containing sample barcodes including the barcode sequence and the sample id. |
| -o OUTPUT, --output OUTPUT | Output file path to populate with the counts for each paired guide and sample. If not provided the counts will be output in STDOUT. |
| --trim-strategy TRIM_STRATEGY | The trim strategy used to extract guides and barcodes. A custom trim strategy should be formatted as as comma separate list of trim coordinates for gRNA1, gRNA2, and the barcode. Each trim coordinate should contain three zero indexed integers giving the file index relative to the order provided in --fastq, the inclusive start index of the trim, and the exclusive end index of the trim. The indices within the trim coordinate should be separated by colon. For convenience the options "two-read" and "three-read" map to default values "0:0:20,1:1:21,1:160:166" and "0:0:20,1:1:21,2:0:6" respectively. The two read strategy should have fastqs R1 and I1. The three read strategy should have fastqs R1, I1, and I2. |
| --gRNA1-error GRNA1_ERROR | The number of substituted base pairs to allow in gRNA1. Must be less than 3. Defaults to 1. |
| --gRNA2-error GRNA2_ERROR | The number of substituted base pairs to allow in gRNA2. Must be less than 3. Defaults to 1. |
| --barcode-error BARCODE_ERROR | The number of insertions, deletions, and subsititions of base pairs to allow in the barcodes. Defaults to 1. |

## Usage example

Test data available in `example-data`.

```
pgmap --fastq example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz \
example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz \
--library example-data/pgPEN-library/paralog_pgRNA_annotations.txt \
--barcodes example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/barcode_ref_file_revcomp.txt \
--trim-strategy two-read
```

For an in-depth example walking through choosing the trim strategy and error tolerances see the guide [here](examples/reproducing-pgpen-two-read-counts.md).

## pgPEN Library

`pgmap` is a general purpose tool for counting reads of paired guide screens, however it was originally developed for the [pgPEN library](https://www.addgene.org/pooled-library/berger-human-pgpen/) described in the paper [Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome](https://pubmed.ncbi.nlm.nih.gov/34469736/) by Parrish _et al_.

## Library Annotations

pgmap searches for guides against a library of known paired guides given as the `--library` argument. This should be a tab separated file with columns labeled "id" "seq_1" and "seq_2". "seq_1" and "seq_2" correspond to the DNA guide sequences (not RNA) that pgmap should map to for each paired guide. "id" should be a unique label for each paired guide, with the suggestion to use include the names of the target genes for each guide. The library annotation for pgPEN is provided in `example-data/pgPEN-library/paralog_pgRNA_annotations.txt`.

## Barcodes

Barcodes should be given in a two column tab separated text file with sample names as the first column and the sample barcode sequences in the second column. An example is provided under `example-data/three-read-strategy/HeLa/screen_barcodes.txt`.

## Custom Trim Strategies

Often the raw sequencing data includes regions of DNA outside of the guides. Trimming these regions to match the expected location of the guides sequences is a critical step to counting guides.

While pgmap provides defaults for a two read and three read trimming strategy to use with pgPEN library screens, it is possible to configure a custom trim strategy using pgmap. A trim strategy is defined by a list of three trim coordinates corresponding to the gRNA1, gRNA2, and barcode locations separated by commas. Each trim coordinate provide a file index zero indexing the location in the list provided in `--fastq`, an inclusive start index for the trim, and an exclusive end index for the trim. Each index should be separated by colons.

For example if two fastqs were provided. The coordinate `1:1:4` would correspond to the second file as given in the order of `--fastq` and would trim to give the second through fourth base pair of the reads.

Here is a full example of the effects of trimming with the default `two-read` trim strategy of `0:0:20,1:1:21,1:160:166`.

Line from input data:

| R1.fastq | I1.fastq |
| - | - |
| **GAGGGAGGGCGAGCTTACGG**TTTTAG | A**CCTCTGGAAGGACACTTCTG**GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGGGTTTGGAGCGAGATTGATAAAGTGCCGACGGTC**CCGGTG**AT |

Corresponding trim output:

| gRNA1 | gRNA2 | barcode |
| - | - | - |
| **GAGGGAGGGCGAGCTTACGG** | **CCTCTGGAAGGACACTTCTG** | **CCGGTG** |

## Output Format

`pgmap` outputs in a tab separated table with the columns for "id", "seq_1", "seq_2" corresponding to the library annotation file and an additional column for counts for each paired guide and sample name defined in the barcodes file.

As an example here is the head of the output from the [usage example](#usage-example) (Note that the counts are zero for these selected guides as the example fastqs are subsetted reads used for functional validation).

| id | seq_1 | seq_2 | sample1 | sample2 | sample3 | sample4 | sample5 |
| - | - | - | - | - | - | - | - |
| AADAC_AADACL2_pg1	| AAGTCTGAAGCACTAAGAAG | AAAGAAAGTCAGAAACCCGA | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg10 | ATTTCTATCCAAATCACTCA | GAAAAAATTTGACTGCAGCA | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg11 | ATTTCTATCCAAATCACTCA | GTGATGTATTCATCTGAAAG | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg12 | ATTTCTATCCAAATCACTCA | TGGGGGCAATTTAGCAACAG | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg13 | GGTATTTCTGGAGATAGTGC | AAAGAAAGTCAGAAACCCGA | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg14 | GGTATTTCTGGAGATAGTGC | GAAAAAATTTGACTGCAGCA | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg15 | GGTATTTCTGGAGATAGTGC | GTGATGTATTCATCTGAAAG | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg16 | GGTATTTCTGGAGATAGTGC | TGGGGGCAATTTAGCAACAG | 0 | 0 | 0 | 0 | 0 |
| AADAC_AADACL2_pg2 | AAGTCTGAAGCACTAAGAAG | GAAAAAATTTGACTGCAGCA | 0 | 0 | 0 | 0 | 0 |
