#' Annotate gimap data
#' @description This is a function here's where we describe what it does
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param annotation_file A special file that contains the list that says which genes should be considered controls or not #TODO: Figure specifics of how this file should be formatted
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#'
#' # Highly recommended but not required
#' run_qc(gimap_dataset)
#'
#' gimap_dataset <- gimap_dataset %>%
#'   gimap_filter() %>%
#'   gimap_annotate()
#'
#' # To see anotations
#' gimap_dataset$annotations
#' }
gimap_annotate <- function(.data = NULL,
                           gimap_dataset,
                           annotation_file = NULL) {

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  if (!is.null(annotation_file)) {
    if (!file.exists(annotation_file)) stop("The annotation_file specified cannot be found. Please double check the file path")
  }

  # TODO: Put the code that annotates the data here!
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/02-get_pgRNA_annotations.Rmd

  gimap_dataset$annotation <- NULL #TODO: Final step is annotations that line up to the same order as the pg gene data should be stored here.

  return(gimap_dataset)
}
