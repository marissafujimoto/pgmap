#' Getting Example Data
#' @description This function pulls some example pgRNA CRIPSR screen data with which to try the tutorials
#' @importFrom readr read_tsv
#' @export
example_data <- function() {
  file <- list.files(
    pattern = "PP_pgPEN_HeLa_counts.txt",
    recursive = TRUE,
    system.file("extdata", package = "metricminer"),
    full.names = TRUE
  )
  readr::read_tsv(file)
}
