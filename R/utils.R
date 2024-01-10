#' Get file path to an key encryption RDS
example_data <- function() {
  file <- list.files(
    pattern = "PP_pgPEN_HeLa_counts.txt",
    recursive = TRUE,
    system.file("extdata", package = "metricminer"),
    full.names = TRUE
  )
  readr::read_tsv(file)
}
