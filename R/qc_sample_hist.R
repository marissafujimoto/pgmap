#' Create a histogram for the pgRNA cpm's, faceted by sample
#' @description This function uses pivot_longer to rearrange the data for plotting and then plots sample specific histograms of the pgRNA cpm's
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @param wide_ar aspect ratio, default is 0.75
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @return sample_cpm_histogram a ggplot
#' @examples \dontrun{
#' 
#' }
#'

qc_sample_hist <- function(gimap_dataset, wide_ar = 0.75){
  
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
  
  return(sample_cpm_histogram)
}