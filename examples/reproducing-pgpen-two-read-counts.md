# Overview

In this document we will demonstrate how `pgmap` is used to produce a table of paired guide counts for a [pgPEN library](https://www.addgene.org/pooled-library/berger-human-pgpen/) CRISPR screen. We demonstrate the ability to use a custom trim strategy to adapt to the sequencing strategy used within the sequencing run. Finally we show how different error tolerances can be configured to maximize the recovery of paired guide counts.

## Data Overview

For this demonstration we have included fastq files, barcodes, and library annotations in `example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5`. More information about the format of these files can be found in the [README](../README.md).

The sequencing run for this screen produced two fastq files: `Undetermined_S0_R1_001_Sampled10k.fastq.gz` and `Undetermined_S0_I1_001_Sampled10k.fastq.gz`. For the purpose of this example, these files have been downsampled to ten thousand reads, but a real sequencing run can have billions of reads.

We also need a barcode file which tells `pgmap` what barcodes correspond to which samples. For this example we provide the barcodes `barcode_ref_file_revcomp.txt` which map to five generic samples: `sample1`, `sample2`, ... `sample5`. In a real experiment the samples may correpond to different time points, treatments, cell lines, or any kind of variable.

### Sequencing Strategy

The sequencing strategy defines how sequences are read from an amplicon. In the case of paired guide CRISPR screens the amplicon contains the CRISPR guide RNAs used for each combinatorial knockout. By reading these guides and counting their abundance or depletion, we can measure the genetic effect of a combinatorial knockout.

### Diagram of Amplicon

![Diagram of pgPEN Library Amplicon: Two Read Sequencing Strategy](assets/two-read-amplicon.svg)

Here we provide a diagram of the amplicon used for pgPEN and a two read strategy.

The R1 read corresponds to the file `Undetermined_S0_R1_001_Sampled10k.fastq.gz`. Similarly the I1 read corresponds to the file `Undetermined_S0_I1_001_Sampled10k.fastq.gz`. Notice how the lengths and start and positions positions of the gRNAs and barcodes differ from the read coordinates. Also note how one read can contain both a guide RNA and a barcode. Our goal is to extract the guide RNAs and barcodes from these reads.

`pgmap` is configurable enough to be used with an arbitrary sequencing strategy. This is accomplished through trim strategies.

## Configuring a Trim Strategy

Using our diagram of the amplicon, we can formulate a trim strategy which extracts the expected locations of each guide RNA and barcode. We say such a location is an _expected_ location since there may be deletions or insertions which offset the guide or barcodes actual position. Alternatively, the sequence could be from controls of the sequencing run (see [PhiX Control v3](https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html)).

When we are configuring `pgmap` with our trim strategy we want to define where we _expect_ the gRNA1, gRNA2, and barcode to be. For the gRNA1, gRNA2, and barcode we define a separate trim coordinate. Each trim coordinate is composed of three integers. The first is a zero-indexed file index, the second a zero-indexed inclusive start index of the read, and finally the last is a zero-indexed exclusive end index of the read.

Using our diagram we can see that gRNA1 is expected to be the first 20 base pairs of the R1 read or using the index format "0:20". Similarly we see that gRNA2 is expected to be the second through 21st read of the I1 read or "1:21". Finally, the barcode is the 161st through 166th base pairs of the I1 read or "160:166".

With our fastq files listed as `Undetermined_S0_R1_001_Sampled10k.fastq.gz`, `Undetermined_S0_I1_001_Sampled10k.fastq.gz` we can point to the R1 file using a 0 and the I1 file using a 1.

Putting this all together we have the trim strategy of "0:0:20,1:1:21,1:160:166". For convenience with the pgPEN library, the trim strategy "two-read" maps to this exact trim strategy.

For more details about the trim strategies and examples of a trim in action see the [README](../README.md).

## Configuring Error Tolerances

When counting guide RNAs and barcodes we want to be able to tolerate errors. These errors could be present before the sequencing (as in a mutant gRNA) or post sequencing (a sequencing read error). `pgmap` can be configured with error tolerances for both gRNA1, gRNA2, and the barcode.

However, there are some nuances with how these error tolerances are counted.

### Guide RNA Error Tolerances

For the guide RNAs the error tolerance is defined using [Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance) or the number of base pair substitutions needed to convert a guide RNA containing an error into the expected guide RNA. There are separate parameters for each guide RNA which can be configured in `pgmap`'s command line interface.

One nuance for the guide RNA error tolerance is that `pgmap` generates every possible misaligned guide RNA in order to efficiently process billions of reads. This means that the higher the error tolerance, the more misaligned guide RNAs need to be generated. Currently `pgmap` supports up to 2 errors in each guide RNA.

### Barcode Error Tolerance

For the barcode, the error tolerance is not defined with Hamming distance, but rather [Levenshtein Distance](https://en.wikipedia.org/wiki/Levenshtein_distance) or the number of base pair deletions, insertions, or substitions to tolerate in the expected barcode location. Using Levenshtein distance allows us to tolerate for insertions and deletions which occur upstream of the barcode trim location.

### Choosing Error Tolerances

In our example, we believe that errors can occur in the sequencing run and we want to be able to account for them. However we don't want to introduce bias in our analysis with the way we are tolerating errors. For our guide RNAs let's keep tolerating a single base pair substitution in either guide since we believe this can happen due to the sequencing technology used (in practice this may be informed by the use of a tool such as [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)). However we want to be careful of choosing the error tolerance for our barcode.

Looking closely at our barcodes we can see that `sample1` has the barcode sequence `CTTGTA` and `sample5` has the sequence `ACTTGA`. Additionally we know that the first base pair of our P7 adaptor was an A. This means that if we expect that deletions could occur in the gRNA2 backbone (see the diagram above), we would read the barcode `CTTGAA` for `sample5` which is only one substitution away from the barcode for `sample1`. Therefore if we don't want to count the cases where there are deletions in the gRNA2 backbone for `sample5` we should not tolerate any errors in barcodes. This may mean that we lose some counts which were legitimately due to errors, but if we are happy with the size of the final counts we can safely ignore these errored reads. Therefore we should use a barcode error tolerance of 0.

The actual error tolerances used in each experiment will depend on the tradeoffs between tolerating errors or strictly focusing on the highest quality sequences only.

## Running `pgmap`

Now we have our files including the fastqs, barcodes, and library annotations. To install `pgmap` refer to the [README](../README.md). In this document we discussed the meaning and choices behind each configuration parameter.

Putting this together in our command line interface we can run the following command:

```
pgmap --fastq example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_R1_001_Sampled10k.fastq.gz \
example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/Undetermined_S0_I1_001_Sampled10k.fastq.gz \
--library example-data/pgPEN-library/paralog_pgRNA_annotations.txt \
--barcodes example-data/two-read-strategy/240123_VH01189_224_AAFGFNYM5/barcode_ref_file_revcomp.txt \
--trim-strategy 0:0:20,1:1:21,1:160:166 \
--gRNA1-error 1 \
--gRNA2-error 1 \
--barcode-error 0 \
```

This should quickly output the counts per paired guide and sample in the format defined in the [README](../README.md). To save this as a file you can either pipe the output directly or use the flag `--output` and pass in a path.

### Notes on Performance

Remember that this example is using a small subsampling of real data and that in practice the runtime of `pgmap` will be much longer. Expect anywhere from ten minutes to a few hours at the extreme higher end to get counts from a large sequencing run. The time taken to run pgmap is roughly proportional to the number of reads in the fastq inputs. Meanwhile the memory used will be proportional to the size of the library in terms of the number of guides and the level of error tolerance desired. `pgmap` was designed to be used in situations where the number of unique guides is on the order of tens of thousands. The expected memory scales linearly with the amount of guides and exponentially with the level of error tolerance desired (max of 2). `pgmap` is single threaded and does not benefit from multithreaded compute resources.
