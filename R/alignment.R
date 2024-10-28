#' Align fastq files to create bam files
#' @description This function takes a forward and reverse bam files for a group of samples and returns the counts and stats.
#' @param fastq_dir a file path to a directory where the bam files associated with the sample names will be stored.
#' @param for_index A file path to the index that was made with `bowtie2_build()` for the forward reference file for the paired guide constructs
#' @param rev_index A file path to the index that was made with `bowtie2_build()` for the reverse reference file for the paired guide constructs
#' @param output_dir A file path that designates where you'd like the aligned files saved. By default this file path is "aligned_bams"
#' @param sample_name a character vector that indicates the name of the samples as they are listed in the file names.
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @param save_logs TRUE/FALSE you want the logs? They will be saved as the sample name
#' @param overwrite TRUE/FALSE should existing bam files of the same name be overwritten?
#' @importFrom furrr future_pmap future_map8p
#' @importFrom purrr reduce
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
#' @examples \dontrun{
#'
#' #################### Build the indices first #################
#'
#' # Declare an input file path for FORWARD
#' for_ref <- file.path(config_folder(), "ref", "paralog_pgRNA1.fa")
#' # Declare ann output file path for FORWARD
#' for_index_file_path <- file.path(config_folder(), "pgPEN_index_for")
#' # Declare an input file path for REVERSE
#' rev_ref <- file.path(config_folder(), "ref", "paralog_pgRNA2.fa")
#' # Declare a output file path for REVERSE
#' rev_index_file_path <- file.path(config_folder(), "pgPEN_index_rev")
#'
#' # Build the FORWARD reference
#' bowtie2_build(
#'   references = for_ref,
#'   bt2Index = for_index_file_path,
#'   overwrite = TRUE
#' )
#'
#' # Build the REVERSE reference
#' bowtie2_build(
#'   references = rev_ref,
#'   bt2Index = rev_index_file_path,
#'   overwrite = TRUE
#' )
#'
#' # This will be the directory that contains all the file path
#' fastq_dir <- file.path(example_data_folder(), "fastqs", "fastq_demuxed")
#'
#' # These MUST be names that are listed in the bam file name itself.
#' # There needs to be exactly 2 bam files per sample
#' sample_names <- c("sample1", "sample2", "sample3")
#'
#' fastq_to_bam(
#'   fastq_dir = fastq_dir,
#'   for_index = for_index_file_path,
#'   rev_index = rev_index_file_path,
#'   sample_names = sample_names,
#'   output_dir = file.path(example_data_folder(), "aligned_bam"),
#'   time = TRUE,
#'   overwrite = TRUE
#' )
#' }
fastq_to_bam <- function(fastq_dir,
                         for_index,
                         rev_index,
                         sample_names,
                         output_dir = "aligned_bam",
                         save_logs = FALSE,
                         time = TRUE,
                         overwrite = TRUE) {
  if (time) timing <- Sys.time()

  # Grab the file names
  sample_df <- grab_paired_files(
    dir = fastq_dir,
    sample_names = sample_names
  ) %>%
    tidyr::pivot_longer(
      cols = c(reverse_file, forward_file),
      names_to = "which_file",
      values_to = "fastq"
    )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Now use the custom function to get the counts for each sample from the pair of bam files
  align <- furrr::future_pmap(sample_df, function(sample_name, fastq, which_file) {
    message("Aligning: ", fastq)

    out_sam_file <- file.path(output_dir, paste0(sample_name, "_", which_file))

    if (which_file == "forward_file") index <- for_index
    if (which_file == "reverse_file") index <- rev_index

    (cmdout <- bowtie2_samtools(
      bt2Index = index,
      output = out_sam_file,
      outputType = "bam",
      seq1 = fastq,
      seq2 = NULL,
      overwrite = overwrite,
      bamFile = NULL
    ))

    out_sam_file <- paste0(out_sam_file, ".bam")

    if (save_logs) writeLines(cmdout, paste0(sample_name, "align_log.txt"))

    return(out_sam_file)
  })

  # Bring along sample names
  names(counts) <- sample_names

  message("Bam files saved to:", output_dir)
}
