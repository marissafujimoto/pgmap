#' Align fastq files
#' @description This function takes a forward and reverse bam files for a group of samples and returns the counts and stats.
#' @param bam_dir a file path to a directory where the bam files associated with the sample names will be stored.
#' @param sample_name a character vector that indicates the name of the samples as they are listed in the file names.
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @importFrom furrr future_pmap future_map
#' @importFrom purrr reduce
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
#' @examples \dontrun{
#'
#' # This will be the directory that contains all the file path
#' fastq_dir <- file.path(example_data_folder(), "fastqs")
#'
#' # These MUST be names that are listed in the bam file name itself.
#' # There needs to be exactly 2 bam files per sample
#' sample_names <- c("sample1", "sample2", "sample3")
#'
#' fastq_to_sam(
#'   fastq_dir = fastq_dir,
#'   sample_names = sample_names,
#'   time = TRUE
#' )
#'
#' }
fastq_to_sam <- function(fastq_dir, index = ,sample_names, time = TRUE) {
  if (time) timing <- Sys.time()
  
  # Grab the file names
  sample_df <- grab_paired_files(
    bam_dir = bam_dir,
    sample_names = sample_names
  )
  
  # Now use the custom function to get the counts for each sample from the pair of bam files
  counts <- furrr::future_pmap(sample_df, function(sample_name, forward_file, reverse_file) {
    sample_count(
      bam_1 = forward_file,
      bam_2 = reverse_file,
      sample_name = sample_name
    )
  })
  # Bring along sample names
  names(counts) <- sample_names
  
  cmdout<-bowtie2_samtools(bt2Index = file.path(example_data_folder(), "pgPEN_index"),
                           output = file.path(td, "result"),
                           outputType = "bam",
                           seq1 = forward_file,
                           seq2 = reverse_file,
                           overwrite = TRUE,
                            "--threads 3"))

  head(readLines(file.path(td, "result.sam")))
}
}
