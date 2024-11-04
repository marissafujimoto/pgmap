#' Obtain the stats for a group of samples
#' @description This function takes a forward and reverse bam files for a group of samples and returns alignment stats.
#' @param bam_dir a file path to a directory where the bam files associated with the sample names will be stored.
#' @param sample_names a character vector that indicates the name of the samples as they are listed in the file names.
#' @param output_dir a file path to where you'd like the stats to be saved to. By default its "stats" folder. This folder will be
#' created if it doesn't exist.
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @importFrom furrr future_map
#' @importFrom Rsamtools quickBamFlagSummary
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
#' stats_df <- get_stats(
#'   bam_dir = bam_dir,
#'   sample_names = sample_names,
#'   output_dir = file.path(example_data_folder(), "bam_stats")
#' )
#'
#' # We can write these to CSV files like this
#' }
get_stats <- function(bam_dir, sample_names, output_dir = "stats", time = FALSE) {
  if (time) timing <- Sys.time()

  # Grab the file names
  sample_df <- grab_paired_files(
    dir = bam_dir,
    sample_names = sample_names
  )

  # Make output directory if it doesn't exist
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

  # Now use the custom function to get the counts for each sample from the pair of bam files
  stats <- furrr::future_pmap(sample_df, function(sample_name, forward_file, reverse_file) {
    # Get stats and save to files
    sink(file.path(output_dir, paste0(sample_name, "_F_stats.txt")), type = "output")
    quickBamFlagSummary(forward_file)
    sink()
    sink(file.path(output_dir, paste0(sample_name, "_F_stats.txt")), type = "output")
    quickBamFlagSummary(forward_file)
    sink()
  })

  if (time) message(paste0("Done(", signif(as.numeric(difftime(Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))

  # Return as a list of two data frames
  message("Stats saved to", output_dir)
}
