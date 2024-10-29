utils::globalVariables(c(
  "X1", "X2", "rname_1", "rname_2", "qname", "paired", "any_paired", "weight",
  "rname", ".", "bam_dir", "forward_file", "reverse_file")
  )

#' Returns file paths to example data for package
#' @description This function loads and returns file paths to example data for the packagae. Which dataset is returned must be specified
#' @param which_data options are "bam" or "fastq"; specifies which example dataset should be returned
#' @export
#' @examples \dontrun{
#'
#' bam_files <- example_data("bam")
#' fastq_files <- example_data("fastq")
#' reference_files <- example_data("ref")
#' }
example_data <- function(which_data) {
  if (which_data == "bam") {
    file <- list.files(
      pattern = ".bam",
      recursive = TRUE,
      system.file("extdata", package = "pgmap"),
      full.names = TRUE
    )
    return(file)
  } else if (which_data == "fastq") {
    file <- list.files(
      pattern = "fastq.gz",
      recursive = TRUE,
      system.file("extdata", package = "pgmap"),
      full.names = TRUE
    )
    return(file)
  } else if (which_data == "ref") {
    file <- list.files(
      pattern = ".fa",
      recursive = TRUE,
      system.file("extdata", package = "pgmap"),
      full.names = TRUE
    )
    return(file)
  } else {
    stop("Specification for `which_data` not understood; Need to use 'bam' or 'fastq'")
  }
}


#' Get file path to example data
#' @export
#' @return Returns the file path to folder where the example data is stored
example_data_folder <- function() {
  file <- list.files(
    pattern = "README.md",
    recursive = TRUE,
    system.file("extdata", package = "pgmap"),
    full.names = TRUE
  )
  file <- grep("idemp", file, invert = TRUE, value = TRUE)
  return(dirname(file))
}

#' Get file path to config folder
#' @export
#' @return Returns the file path to folder where the example data is stored
config_folder <- function() {
  file.path(example_data_folder(), "..", "config")
}
