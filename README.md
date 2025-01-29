# pgmap - (paired guide RNA MAPper) Pipeline

We developed pgMAP, an analysis pipeline to map gRNA sequencing reads from dual-targeting CRISPR screens. pgMAP output includes a dual gRNA read counts table and quality control metrics including the proportion of correctly-paired reads and CRISPR library sequencing coverage across all time points and samples.

It is based off of the original code and research from the Berger Lab stored in this repository: [https://github.com/FredHutch/GI_mapping](https://github.com/FredHutch/pgMAP_pipeline)

## Prerequisites

pgmap requires python3 above or equal to version 3.10

### Installation using pip

```
pip install pgmap
```

Check if pgmap is installed in path:
```
pgmap --help
```

### Installation from source code

First you can clone this repo:
```
git clone https://github.com/FredHutch/pgmap
```

To install required dependencies:
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
- barcodes path which says which barcodes go to which cells
- pgPEN annotations path (this is included in the package)

Usage example (test data available in `example-data`)
```
pgmap --fastq example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz --library example-data/pgPEN-library/paralog_pgRNA_annotations.txt --barcodes example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/barcode_ref_file_revcomp.txt --trim-strategy two-read
```
