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

# TODO: Candace looks into what kind of R class we want
gimap_data <- R6::R6Class("gimap_data",
            public = list(
              metadata = NULL,
              counts = NULL,
              counts_per_sample = NULL
            ))

#' Making a new gimap dataset
#' @description This function allows people to have their data ready to be processed by gimap
#' @param metadata
#' @param counts
#' @export
#' @examples \dontrun{
#' data <- example_data() %>%
#'   dplyr::select(c("Day00_RepA", "Day05_RepA", "Day22_RepA", "Day22_RepB", "Day22_RepC")) %>%
#'   as.matrix()
#'
#' counts_data <- setup_data(data)
#' }
setup_data <- function(metadata = NULL, counts = NULL) {

  new_data <- gimap_data$new()

  if (!is.null(counts)) stop("counts cannot be NULL")
  if (!is.matrix(counts)) stop("counts can only be in the form of a matrix")

  new_data$counts <- counts

  new_data$counts_per_sample <- apply(counts, 2, sum)

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) stop("metadata can only be in the form of a data.frame")

    if (!nrow(metadata) != ncol(counts)) {
      stop("the number of rows in the metadata is not equal to the number of columns in counts.
           Are you sure this metadata goes with these counts?")
    }
    new_data$metadata <- metadata
  }

  return(new_data)
}

