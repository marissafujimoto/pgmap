#'
#'
#'
#'
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @param cutoff default is NULL
#' @param wide_ar aspect ratio, default is 0.75
#' @importFrom magrittr %>%
#' @import ggplot2
#' Assumes the first column of the dataset is the day0 data; do I need a better
#' method to tell, especially if there are reps?

qc_plasmid_histogram <- function(gimap_dataset, cutoff = NULL, wide_ar = 0.75){
  to_plot <- data.frame(gimap_dataset$transformed_data$log2_cpm[,1]) %>% `colnames<-`(c("log2_cpm"))
  plasmid_cpm_histogram <- ggplot(to_plot, aes(x=log2_cpm)) +
    geom_histogram(binwidth=0.2, color = "black", fill="gray60") +
    plot_options() +
    plot_theme() +
    theme(aspect.ratio = wide_ar)
  
  if (is.null(cutoff)){
    #if cutoff is null, suggest a cutoff and plot with suggested
    quantile_info <- quantile(to_plot$log2_cpm)
    plasmid_cpm_stats <- data.frame(stat = c("median", "Q1", "Q3", "lower_outlier"),
                                    log2_cpm_value = c(quantile_info["50%"],
                                              quantile_info["25%"],
                                              quantile_info["75%"],
                                              quantile_info["25%"] - (1.5*(quantile_info["75%"] - quantile_info["25%"]))
                                    )
    )
    cutoff <- plasmid_cpm_stats[which(plasmid_cpm_stats$stat == "lower_outlier"),"log2_cpm_value"]
  } else{ plasmid_cpm_stats <- NULL }
  
  #plot with the cutoff
  plasmid_cpm_hist_wcutoff <- plasmid_cpm_histogram +
    geom_vline(xintercept = cutoff,
               linetype = "dashed")
  
  
  plasmid_cpm_filter <- unlist(lapply(1:nrow(to_plot), function(x) to_plot$log2_cpm[x] < cutoff))
  
  plasmid_filter_df <- data.frame("Plasmid_log2cpmBelowCutoff" = c(FALSE, TRUE), n = c(sum(!plasmid_cpm_filter), sum(plasmid_cpm_filter))) %>%
    mutate(percent = round(((n/sum(n))*100),2))
  
  return(list(plasmid_hist_nocutoff= plasmid_cpm_histogram, 
              plasmid_stats <- plasmid_cpm_stats,
              used_log2_cpm_cutoff = cutoff,
              plasmid_hist_cutoff=plasmid_cpm_hist_wcutoff,
              plasmid_filter = plasmid_cpm_filter,
              plasmid_filter_report = plasmid_filter_df))
}