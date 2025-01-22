# pgmap - (paired guide RNA MAPper) Pipeline

We developed pgMAP, an analysis pipeline to map gRNA sequencing reads from dual-targeting CRISPR screens. pgMAP output includes a dual gRNA read counts table and quality control metrics including the proportion of correctly-paired reads and CRISPR library sequencing coverage across all time points and samples.

It is based off of the original code and research from the Berger Lab stored in this repository: [https://github.com/FredHutch/GI_mapping](https://github.com/FredHutch/pgMAP_pipeline)

## Prerequisites

In order to run this pipeline you will need R and to install the `pgmap` package and its dependencies. In R you can run this to install the package:


At this time, pgmap is not yet published to PyPI. But can be used following these steps:

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
python3 -m pip install .
```

## Getting Started

You'll need to declare four files (we have example data to work with):

- Two sequencing reads files:
  - I1 fastq file
  - R1 fastq file
- barcodes path which says which barcodes go to which cells
- pgPEN annotations path (this is included in the package)
