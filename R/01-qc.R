
#' This is a title for a function
#' @description This is a function here's where we describe what it does
#' @param parameter Here's a parameter let's describe it here
#' @export
#' @importFrom tidyr pivot_longer
#' @import ggplot2 pheatmap
#' @examples \dontrun{
#'
#' }


run_qc <- function(gimap_dataset, plots_dir = "./qc_plots", wide_ar = 0.75, square_ar = 1) {
  
  if (!dir.exists(plots_dir)){
    dir.create(plots_dir, showWarnings = TRUE)
  }
  
  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")
  
  long_form <- 
    tidyr::pivot_longer(data.frame(gimap_dataset$transformed_data$count_norm),
                        everything(),
                        names_to = "sample",
                        values_to = "count_normalized")
  
  counts_cdf <- ggplot(long_form, aes(x=count_normalized, color = sample)) +
    stat_ecdf() +
    labs(x = "-log10(count/total_count)", # bquote(~-log[10]~"(count/total_count)")
         y = "Expected_pgRNAs",
         color = "Sample") +
    plot_options() +
    plot_theme() +
    theme(aspect.ratio = wide_ar)
  
  counts_cdf
  save_plot(counts_cdf)
  
  long_form <-
    tidyr::pivot_longer(data.frame(gimap_dataset$transformed_data$log2_cpm),
                        everything(),
                        names_to = "sample",
                        values_to = "log2_cpm")
  
  sample_cpm_histogram <- ggplot(long_form, aes(x=log2_cpm, fill = sample)) +
    geom_histogram(color = "black", binwidth = 0.5) +
    plot_options() +
    plot_theme() +
    theme(aspect.ratio = wide_ar,
          legend.position = "none") +
    facet_wrap(~sample, scales = "free_y", ncol = ceiling(ncol(gimap_dataset$raw_counts)/2))
  
  sample_cpm_histogram
  save_plot(sample_cpm_histogram)
  
  cpm_cor <- gimap_dataset$transformed_data$cpm %>%
    cor() %>%
    round(2) %>%
    data.frame()
  
  sample_cor_heatmap_unfiltered <-
    pheatmap(cpm_cor,
             border_color = "white",
             cellwidth = 20, cellheight = 20,
             treeheight_row = 20, treeheight_col = 20)
  
  sample_cor_heatmap_unfiltered
  
  ## flag pgRNAs with count = 0 at any time point
  counts_cpm_filter <- unlist(lapply(1:nrow(gimap_dataset$raw_counts), function(x) 0 %in% gimap_dataset$raw_counts[x,]))
  
  ## how many pgRNAs have a raw count of 0 at any time point/sample
  d.summ <- data.frame("RawCount0" = c(FALSE, TRUE), n = c(sum(!counts_cpm_filter), sum(counts_cpm_filter))) %>%
    mutate(percent = round(((n/sum(n))*100),2))
  
  d.summ
  
  ## flag pgRNAs with low plasmid read counts
    
  
}
