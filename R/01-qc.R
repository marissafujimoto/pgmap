
#' A function to run QC
#' @description This is a function here's where we describe what it does
#' @param parameter Here's a parameter let's describe it here
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @import ggplot2 pheatmap
#' @examples \dontrun{
#'
#' }


run_qc <- function(gimap_dataset, plasmid_cutoff = NULL, use_combined = TRUE,
                   params = list(cell_line = NULL), plots_dir = "./qc_plots", overwrite = FALSE, output_format = NULL, output_file = "QC_Report", output_dir = "./",
                   wide_ar = 0.75, square_ar = 1) {
  
  if (!dir.exists(plots_dir)){
    dir.create(plots_dir, showWarnings = TRUE)
  }
  
  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")
  
  # for Rmd report
  list_of_qc_things <- list(counts_cdf = NULL,
                            sample_cpm_histogram = NULL,
                            sample_cor_heatmap_unfiltered = NULL,
                            zerocount_df = NULL,
                            plasmid_hist_nocutoff = NULL,
                            plasmid_stats = NULL,
                            plasmid_hist_cutoff = NULL,
                            plasmid_filter_report = NULL,
                            combined_filter_df = NULL,
                            which_filter_df = NULL
                            )
  
  counts_cdf <- qc_cdf(gimap_dataset, wide_ar = wide_ar)
  
  list_of_qc_things$counts_cdf <- counts_cdf #put it in the report
  save_plot(counts_cdf, params = params, out_dir = plots_dir)
  
  sample_cpm_histogram <- qc_sample_hist(gimap_dataset, wide_ar = wide_ar)
  
  list_of_qc_things$sample_cpm_histogram <- sample_cpm_histogram #put it in the report
  save_plot(sample_cpm_histogram, params = params, out_dir = plots_dir)
  
  sample_cor_heatmap_unfiltered <- qc_cor_heatmap(gimap_dataset)
  
  list_of_qc_things$sample_cor_heatmap_unfiltered <- sample_cor_heatmap_unfiltered #put it in the report
  save_plot(sample_cor_heatmap_unfiltered, params = params, out_dir = plots_dir)
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
  
  list_of_qc_things$plasmid_hist_nocutoff <- plasmid_filter_list$plasmid_hist_nocutoff #put it in the report
  save_plot(plasmid_filter_list$plasmid_hist_nocutoff, params = params, out_dir = plots_dir)
  
  
  list_of_qc_things$plasmid_stats <- plasmid_filter_list$plasmid_stats #put it in the report

  list_of_qc_things$plasmid_hist_cutoff <- plasmid_filter_list$plasmid_hist_cutoff #put it in the report
  save_plot(plasmid_filter_list$plasmid_hist_cutoff, params = params, out_dir = plots_dir)
  
  
  list_of_qc_things$plasmid_filter_report <- plasmid_filter_list$plasmid_filter_report #put it in the report
  
  plasmid_filter <- plasmid_filter_list$plasmid_filter #TRUE means the filter applies
  
  ## combine the filters (zerocount and plasmid)
  ## how many pgRNAs are removed by both filters?
  
  combined_filter_list <- qc_combine_filters(counts_filter, plasmid_filter, use_combined = use_combined)
  
  list_of_qc_things$num_filtered_report <- combined_filter_df <- combined_filter_list$num_filtered_report #put it in the report
  
  list_of_qc_things$which_filter_df <- which_filter_df <- combined_filter_list$num_which_filter_report #put it in the report
  
  gimap_dataset$transformed_data$qc_filter <- combined_filter_list$combined_filter
  
  #if (use_combined){
    #filter gimap_dataset and metadata to only keep pgRNA passing both filters
    #gimap_dataset$transformed_data$qc_filter <- combined_filter_list$keep_pgRNA
  #} else {
    #print out a message asking which filter to use?
  #}
  
  qc_generateReport(list_of_qc_things, overwrite = overwrite, output_file = output_file, output_dir = output_dir)
  
  return(gimap_dataset)
}
