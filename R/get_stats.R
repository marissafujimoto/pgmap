#' Obtain the counts for a group of samples
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
#' bam_dir <- file.path(example_data_folder(), "bam") 
#' 
#' # These MUST be names that are listed in the bam file name itself.
#' # There needs to be exactly 2 bam files per sample
#' sample_names <- c("sample1", "sample2", "sample3")
#' 
#' counts_df <- calc_counts(
#'   bam_dir = bam_dir, 
#'   sample_names = sample_names, 
#'   time = TRUE
#')
#'
#' # We can write these to CSV files like this
#' readr::write_csv(counts_df$counts, "counts_pgmap.tsv")
#' readr::write_csv(counts_df$stats, "stats_pgmap.tsv")
#' 
#' }
get_stats <- function(bam_dir, sample_names, time = TRUE) {
  
  if (time) timing <- Sys.time()
  
  # Grab the file names
  sample_df <- grab_paired_files(bam_dir = bam_dir, 
                                 sample_names = sample_names)
  
  # Now use the custom function to get the counts for each sample from the pair of bam files
  stats <- furrr::future_pmap(sample_df, function(sample_name, forward_file, reverse_file) {
    # Read in the BAM files
    bam_1_stats <- quickBamFlagSummary(forward_file)
    bam_2_stats <- quickBamFlagSummary(reverse_file)
  })
  # Bring along sample names
  names(counts) <- sample_names
  
  # Then condense all samples to one stats data frame
  stats_df <- furrr::future_map(counts, "stats") %>%
    dplyr::bind_rows(.id = "sample")
  
  # Then condense all samples to one counts data frame
  counts_df <- furrr::future_map(counts, "counts", right_join) %>%
    purrr::reduce(dplyr::full_join)
  
  if (time) message(paste0("Done(", signif(as.numeric(difftime (Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))
  # Return as a list of two data frames
  return(list(counts = counts_df, stats = stats_df))
}

#' Obtain the counts for a single sample
#' @description This function takes a forward and reverse bam file and returns the counts
#' @param bam_1 a file path to one of the bam files for the sample
#' @param bam_2 a file path to the other one of the bam files for the sample
#' @param sample_name a single character that indicates the name of the sample "sample1"
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @import dplyr
#' @importFrom Rsamtools ScanBamParam scanBam
#' @export
#' @examples \dontrun{
#'
#' bam_dir <- file.path(example_data_folder(), "bam") 
#' 
#' counts <- sample_stats(
#'   bam_1 = file.path(bam_dir, "pgMAP_tutorial_gRNA1_trimmed_sample1_aligned.bam"),
#'   bam_2 = file.path(bam_dir, "pgMAP_tutorial_gRNA2_trimmed_sample1_aligned.bam"),
#'   sample_name = "sample1", 
#'   time = TRUE)
#'
#' }

sample_stats <- function(bam_1, bam_2, sample_name, time = FALSE) {
  
  if (time) timing <- Sys.time()
  # Print out message
  message(paste("Getting stats for", sample_name))
  
  # Read in the BAM files
  bam_1 <- quickBamFlagSummary(bam_1, param = param)
  bam_2 <- quickBamFlagSummary(bam_2, param = param)
  
  if (time) message(paste0("Done(", signif(as.numeric(difftime (Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))
  return(list(counts = counts_df, stats = stats_df))
}

