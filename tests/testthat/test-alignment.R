test_that("Test alignment", {
  
   # Declare an input file path for FORWARD
   for_ref <- file.path(config_folder(), "ref", "paralog_pgRNA1.fa")
   # Declare ann output file path for FORWARD
   for_index_file_path <- file.path(config_folder(), "pgPEN_index_for")
   # Declare an input file path for REVERSE
   rev_ref <- file.path(config_folder(), "ref", "paralog_pgRNA2.fa")
   # Declare a output file path for REVERSE
   rev_index_file_path <- file.path(config_folder(), "pgPEN_index_rev")
   
   # Build the FORWARD reference
   bowtie2_build(
     references = for_ref,
     bt2Index = for_index_file_path,
     overwrite = TRUE
   )
   
   # Build the REVERSE reference
   bowtie2_build(
     references = rev_ref,
     bt2Index = rev_index_file_path,
     overwrite = TRUE
   )
  
   # This will be the directory that contains all the file path
   fastq_dir <- file.path(example_data_folder(), "fastqs", "fastq_demuxed")
  
   # These MUST be names that are listed in the bam file name itself.
   # There needs to be exactly 2 bam files per sample
   sample_names <- c("sample1", "sample2", "sample3")
  
   fastq_to_bam(
     fastq_dir = fastq_dir,
     for_index = for_index_file_path,
     rev_index = rev_index_file_path,
     sample_names = sample_names,
     output_dir = file.path(example_data_folder(), "aligned_bam"),
     time = TRUE, 
     overwrite = TRUE
   )
})
