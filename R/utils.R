#' This is a title for a function
#' @description This is a function here's where we describe what it does
#' @param parameter Here's a parameter let's describe it here
#' @export
#' @examples \dontrun{
#'
#' pg_data <- example_data()
#' }
example_data <- function() {
  file <- list.files(
    pattern = "PP_pgPEN_HeLa_counts.txt",
    recursive = TRUE,
    system.file("extdata", package = "gimap"),
    full.names = TRUE
  )
  readr::read_tsv(file)
}
