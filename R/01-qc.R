#' A function to run QC
#' @description This is a function here's where we describe what it does
#' @param plasmid_cutoff default it NULL; `qc_plasmid_histogram` will calculate and set a default one unless a value is specified here
#' @param use_combined default it TRUE; if TRUE, both zero count and low plasmid CPM filters are applied and if either is TRUE, a pgRNA construct will be filtered out. If FALSE, need to allow to specify which should be used
#' @param plots_dir default is `./qc_plots`; directory to save plots created with this function, if it doesn't exist already it will be created
#' @param overwrite default is FALSE; whether to overwrite the QC Report file
#' @param output_file_path default is `QC_Report`; name of the output QC report file
#' @param ...
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @examples \dontrun{
#'
#' }
run_qc <- function(gimap_dataset,
                   plots_dir = "./qc_plots",
                   overwrite = FALSE,
                   output_format = NULL,
                   output_file = "QC_Report",
                   output_dir = "./",
                   plot_suffix = "png") {

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, showWarnings = TRUE)
  }

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  # Determine the template
  templateFile <- system.file("rmd/gimapQCTemplate.Rmd", package = "gimap")

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
    stop("The Rmd template file ", templateFile, " does not exist", call. = FALSE)
  }

  rmarkdown::render(params = list(dataset = gimap_dataset,
                                  plots_dir = plots_dir),
                                  ...)

  return(gimap_dataset)
}
