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
#' sample_df <- grab_paired_files(bam_dir, sample_names)
#' }
grab_paired_files <- function(bam_dir, sample_names, time = TRUE) {
  # Get the file paths for each sample name
  sample_files <- sapply(
    sample_names, function(file_name) {
      files <- list.files(path = bam_dir, pattern = file_name, full.names = TRUE)

      if (length(files) < 2) stop(file_name, ": does not have 2 files associated with it.")
      if (length(files) > 2) stop(file_name, ": has more than 2 files associated with it.")

      return(files)
    },
    USE.NAMES = TRUE
  )

  # Reformat so samples are rows and files are in columns
  sample_df <- sample_files %>%
    t() %>%
    data.frame() %>%
    dplyr::rename(
      forward_file = X1,
      reverse_file = X2
    ) %>%
    tibble::rownames_to_column("sample_name")

  return(sample_df)
}
