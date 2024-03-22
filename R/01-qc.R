#' A function to run QC
#' @description This is a function here's where we describe what it does
#' @param use_combined default it TRUE; if TRUE, both zero count and low plasmid CPM filters are applied and if either is TRUE, a pgRNA construct will be filtered out. If FALSE, need to allow to specify which should be used
#' @param plots_dir default is `./qc_plots`; directory to save plots created with this function, if it doesn't exist already it will be created
#' @param overwrite default is FALSE; whether to overwrite the QC Report file
#' @param output_file_path default is `QC_Report`; name of the output QC report file
#' @param ... additional parameters are sent to rmarkdown::render
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @examples \dontrun{
#'
#' }
run_qc <- function(gimap_dataset,
                   output_file = "./gimap_QC_Report.Rmd",
                   plots_dir = "./qc_plots",
                   overwrite = FALSE,
                   ...) {

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  # Determine the template
  templateFile <- system.file("rmd/gimapQCTemplate.Rmd", package = "gimap")

  # Make sure that the directory exists!
  directory <- dirname(output_file)

  # If not, make the directory
  if (directory != ".") {
    if (!dir.exists(directory)) dir.create(directory, showWarnings = TRUE, recursive = TRUE)
  }

  # Now if the file exists,
  if (file.exists(templateFile)) {
    if (file.exists(output_file) & !overwrite) {
      stop("there is already an output .Rmd file", output_file,
           ". Please remove or rename this file, or choose another output_file name.",
           call. = FALSE
      )
    } else {
      file.copy(from = templateFile, to = output_file, overwrite = overwrite)
    }
  } else {
    stop("The Rmd template file ", templateFile, " does not exist -- did you move it from the package files?",
         call. = FALSE)
  }

  # Make a plots directory if it doesn't exist
  if (!dir.exists(plots_dir)) dir.create(plots_dir, showWarnings = TRUE)

  # Send the data to render it!
  rmarkdown::render(output_file,
                    params = list(dataset = gimap_dataset,
                                  plots_dir = plots_dir),
                                  ...)

  # Tell where the output is
  results_file <- gsub("\\.Rmd$", "\\.html", output_file)
  message("Results in: ", results_file)

  results_file <- normalizePath(list.files(pattern = results_file, full.names = TRUE))

  if (results_file != "") browseURL(results_file)
}
