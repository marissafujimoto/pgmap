#' Create a CDF for the pgRNA normalized counts
#' @description This function uses pivot_longer to rearrange the data for plotting and then plots a CDF of the normalized counts
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @param wide_ar aspect ratio, default is 0.75
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @return counts_cdf a ggplot
#' @examples \dontrun{
#' 
#' }
#'

qc_cdf <- function(gimap_dataset, wide_ar = 0.75){
  
  long_form <- 
    tidyr::pivot_longer(data.frame(gimap_dataset$transformed_data$count_norm),
                        everything(),
                        names_to = "sample",
                        values_to = "count_normalized")
  
  counts_cdf <- ggplot(long_form, aes(x=count_normalized, color = sample)) +
    stat_ecdf() +
    labs(x = "-log10(count/total_count)",
         y = "Expected_pgRNAs",
         color = "Sample") +
    plot_options() +
    plot_theme() +
    theme(aspect.ratio = wide_ar)
  
  return(counts_cdf)
}