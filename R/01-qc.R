
#' A function to run QC
#' @description This is a function here's where we describe what it does
#' @param plasmid_cutoff default it NULL; `qc_plasmid_histogram` will calculate and set a default one unless a value is specified here
#' @param use_combined default it TRUE; if TRUE, both zero count and low plasmid CPM filters are applied and if either is TRUE, a pgRNA construct will be filtered out. If FALSE, need to allow to specify which should be used
#' @param params default is `list(cell_line = NULL)`; used to pass info when saving plots
#' @param plots_dir default is `./qc_plots`; directory to save plots created with this function, if it doesn't exist already it will be created
#' @param overwrite default is FALSE; whether to overwrite the QC Report file
#' @param output_format default is NULL; the output format that the QC Rmd will be knitted to/produced
#' @param output_file default is `QC_Report`; name of the output QC report file
#' @param output_dir default is `./`; where the qc report will be stored
#' @param wide_ar default is 0.75; aspect ratio
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @examples \dontrun{
#'
#' }


run_qc <- function(gimap_dataset,
                   plasmid_cutoff = NULL,
                   use_combined = TRUE,
                   params = list(cell_line = NULL),
                   plots_dir = "./qc_plots",
                   overwrite = FALSE,
                   output_format = NULL,
                   output_file = "QC_Report",
                   output_dir = "./",
                   wide_ar = 0.75,
                   plot_suffix = "png") {


  if (!dir.exists(plots_dir)){
    dir.create(plots_dir, showWarnings = TRUE)
  }

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  # for Rmd report
  qc_obj <- list(counts_cdf = NULL,
                            counts_cdf_file_loc = NULL,

                            sample_cpm_histogram = NULL,
                            sample_cpm_histogram_file_loc = NULL,

                            sample_cor_heatmap_unfiltered = NULL,
                            sample_cor_heatmap_unfilered_file_loc = NULL,

                            zerocount_df = NULL,
                            plasmid_hist_nocutoff = NULL,
                            plasmid_hist_nocutoff_file_loc = NULL,

                            plasmid_stats = NULL,
                            plasmid_hist_cutoff = NULL,
                            plasmid_filter_report = NULL,
                            num_filtered_report = NULL,
                            which_filter_df = NULL
                            )

  qc_obj <- qc_cdf(gimap_dataset,
                   qc_obj,
                   plots_dir = plots_dir,
                   wide_ar = wide_ar)

  qc_obj <- qc_sample_hist(gimap_dataset,
                           qc_obj,
                           plots_dir = plots_dir,
                           wide_ar = wide_ar)

  qc_obj <- qc_cor_heatmap(gimap_dataset,
                           qc_obj,
                           plots_dir = plots_dir)

  qc_obj <- qc_filter_zerocounts(gimap_dataset,
                                 qc_obj,
                                 plots_dir = plots_dir)

  qc_obj <- qc_plasmid_histogram(gimap_dataset,
                                 qc_obj,
                                 plots_dir = plots_dir,
                                 cutoff = plasmid_cutoff,
                                 wide_ar = wide_ar)

  qc_obj <- qc_combine_filters(counts_filter,
                               plasmid_filter,
                               use_combined = use_combined)

  qc_generateReport(qc_obj,
                    overwrite = overwrite,
                    output_file = output_file,
                    output_dir = output_dir)

  return(gimap_dataset)
}
