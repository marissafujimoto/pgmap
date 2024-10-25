#' Align fastq files to create bam files
#' @description This function takes a forward and reverse bam files for a group of samples and returns the counts and stats.
#' @param bam_dir a file path to a directory where the bam files associated with the sample names will be stored.
#' @param sample_name a character vector that indicates the name of the samples as they are listed in the file names.
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @param save_logs TRUE/FALSE you want the logs? They will be saved as the sample name
#' @importFrom furrr future_pmap future_map
#' @importFrom purrr reduce
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
#' @examples \dontrun{
#'
#' # Build the index first
#' refs <- file.path(config_folder(), "ref", c("pgPEN_R1.fa", "pgPEN_R2.fa"))
#'
#' cmdout <- bowtie2_build(
#'   references = refs,
#'   bt2Index = file.path(config_folder(), "pgPEN_index"),
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
#'   sample_names = sample_names,
#'   output_dir = file.path(example_data_folder(), "sam"),
#'   time = TRUE, 
#'   overwrite = TRUE
#' )
#' }
fastq_to_bam <- function(fastq_dir,
                         index = file.path(config_folder(), "pgPEN_index"),
                         sample_names,
                         sort = TRUE,
                         output_dir = "sam",
                         save_logs = FALSE,
                         time = TRUE, 
                         overwrite = TRUE) {
  if (time) timing <- Sys.time()

  # Grab the file names
  sample_df <- grab_paired_files(
    dir = fastq_dir,
    sample_names = sample_names
  ) %>% 
   tidyr::pivot_longer(cols = c(reverse_file, forward_file), 
                       names_to = "which_file", 
                       values_to = "fastq")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Now use the custom function to get the counts for each sample from the pair of bam files
  align <- furrr::future_pmap(sample_df, function(sample_name, fastq, which_file) {
    
    message("Aligning:", sample_name, " ", which_file)
    
    out_sam_file <- file.path(output_dir, paste0(sample_name, "_", which_file))
                          
    (cmdout <- bowtie2_samtools(
      bt2Index = index,
      output = out_sam_file,
      outputType = "sam",
      seq1 = fastq,
      seq2 = NULL,
      overwrite = overwrite,
      bamFile = NULL, 
      "--quiet"
    ))
    
    out_sam_file <- paste0(out_sam_file, ".sam")
    
    message("Converting to BAM")
    
    asBam(file = out_sam_file, 
          files = gsub(".sam$", ".bam", out_bam_file), 
          overwrite = overwrite)

    if (save_logs) writeLines(cmdout, paste0(sample_name, "align_log.txt"))

    return(cmdout)
  })

  # Bring along sample names
  names(counts) <- sample_names

  message("Bam files saved to ")
}
