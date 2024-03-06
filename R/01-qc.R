
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


run_qc <- function(gimap_dataset, plasmid_cutoff = NULL, use_combined = TRUE,
                   params = list(cell_line = NULL), plots_dir = "./qc_plots", overwrite = FALSE, output_format = NULL, output_file = "QC_Report", output_dir = "./",
                   wide_ar = 0.75) {
  
  store_results <- function(result_obj, lvoi, list_of_qc_things, plots_dir, params){
    pdf_base <- file.path(plots_dir, "plots", "pdf")
    png_base <- file.path(plots_dir, "plots", "png")
    list_of_qc_things[[lvoi]] <- result_obj
    save_plot(result_obj, params=params, out_dir = plots_dir)
    list_of_qc_things[[paste0(eval(lvoi), "_file_loc")]][["pdf"]] <- file.path(pdf_base, paste0(params$cell_line, "_", deparse(substitute(result_obj)), ".pdf"))
    list_of_qc_things[[paste0(eval(lvoi), "_file_loc")]][["png"]] <- file.path(png_base, paste0(params$cell_line, "_", deparse(substitute(result_obj)), ".png"))
    return(list_of_qc_things)
  }
  
  if (!dir.exists(plots_dir)){
    dir.create(plots_dir, showWarnings = TRUE)
  }
  
  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")
  
  # for Rmd report
  list_of_qc_things <- list(counts_cdf = NULL,
                            counts_cdf_file_loc = list(pdf = NULL,
                                                       png = NULL),
                            sample_cpm_histogram = NULL,
                            sample_cpm_histogram_file_loc = list(pdf = NULL,
                                                                 png = NULL),
                            sample_cor_heatmap_unfiltered = NULL,
                            sample_cor_heatmap_unfilered_file_loc = list(pdf = NULL,
                                                                         png = NULL),
                            zerocount_df = NULL,
                            plasmid_hist_nocutoff = NULL,
                            plasmid_hist_nocutoff_file_loc = list(pdf = NULL,
                                                                  png = NULL),
                            plasmid_stats = NULL,
                            plasmid_hist_cutoff = NULL,
                            plasmid_filter_report = NULL,
                            num_filtered_report = NULL,
                            which_filter_df = NULL
                            )
  
  counts_cdf <- qc_cdf(gimap_dataset, wide_ar = wide_ar)
  
  list_of_qc_things <- store_results(counts_cdf, "counts_cdf", list_of_qc_things, plots_dir, params) #put it in the report
  
  sample_cpm_histogram <- qc_sample_hist(gimap_dataset, wide_ar = wide_ar)
  
  list_of_qc_things <- store_results(sample_cpm_histogram, "sample_cpm_histogram", list_of_qc_things, plots_dir, params) #put it in the report
  
  sample_cor_heatmap_unfiltered <- qc_cor_heatmap(gimap_dataset)
  
  list_of_qc_things <- store_results(sample_cor_heatmap_unfiltered, "sample_cor_heatmap_unfiltered", list_of_qc_things, plots_dir, params) #put it in the report
  ## repeat after filtering?
  
  ## flag pgRNAs with count = 0 at any time point
  ## how many pgRNAs have a raw count of 0 at any time point/sample
  counts_filter_list <- qc_filter_zerocounts(gimap_dataset)
  
  counts_filter <- counts_filter_list$filter #TRUE means the filter applies
  
  zerocount_df <- counts_filter_list$reportdf
  
  list_of_qc_things$zerocount_df <- zerocount_df #put it in the report
  
  ## flag pgRNAs with low plasmid read counts
  ## plasmid reads are the first sample/day0 time point  
  plasmid_filter_list <- qc_plasmid_histogram(gimap_dataset, cutoff = plasmid_cutoff, wide_ar = wide_ar)
  list_of_qc_things$plasmid_stats <- plasmid_filter_list$plasmid_stats #put it in the report
  list_of_qc_things$plasmid_filter_report <- plasmid_filter_list$plasmid_filter_report #put it in the report
  plasmid_filter <- plasmid_filter_list$plasmid_filter #TRUE means the filter applies
  
  plasmid_hist_nocutoff <- plasmid_filter_list$plasmid_hist_nocutoff
  list_of_qc_things <- store_results(plasmid_hist_nocutoff, "plasmid_hist_nocutoff", list_of_qc_things, plots_dir, params) #put it in the report
  
  plasmid_hist_cutoff <- plasmid_filter_list$plasmid_hist_cutoff
  list_of_qc_things <- store_results(plasmid_hist_cutoff, "plasmid_hist_cutoff", list_of_qc_things, plots_dir, params) #put it in the report
  
  ## combine the filters (zerocount and plasmid)
  ## how many pgRNAs are removed by both filters?
  
  combined_filter_list <- qc_combine_filters(counts_filter, plasmid_filter, use_combined = use_combined)
  
  list_of_qc_things$num_filtered_report <- combined_filter_list$num_filtered_report #put it in the report
  
  list_of_qc_things$which_filter_df <- combined_filter_list$num_which_filter_report #put it in the report
  
  gimap_dataset$qc$qc_filter <- combined_filter_list$combined_filter
  
  #if (use_combined){
    #filter gimap_dataset and metadata to only keep pgRNA passing both filters
    #gimap_dataset$qc$qc_filter <- combined_filter_list$keep_pgRNA
  #} else {
    #print out a message asking which filter to use?
  #}
  
  qc_generateReport(list_of_qc_things, overwrite = overwrite, output_file = output_file, output_dir = output_dir)
  
  return(gimap_dataset)
}
