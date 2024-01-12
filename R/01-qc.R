
#' This is a title for a function
#' @description This is a function here's where we describe what it does
#' @param parameter Here's a parameter let's describe it here
#' @export
#' @importFrom tidyr pivot_longer
#' @import ggplot
#' @examples \dontrun{
#'
#' }

run_qc <- function(gimap_data, plots_dir = "./qc_plots", wide_ar = 0.75, square_ar = 1) {

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, showWarnings = TRUE)
  }

  if (!("gimap_dataset" %in% class(gimap_data))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  long_form <-
    tidyr::pivot_longer(data.frame(gimap_dataset$transformed_data$count_norm),
                        everything(),
                        names_to = "sample",
                        values_to = "count normalized")

  counts_cdf <- ggplot(long_form, aes(x = count, color = sample)) +
    stat_ecdf() +
    labs(x = "-log10(count/total_count)", # bquote(~-log[10]~"(count/total_count)")
         y = "Expected_pgRNAs",
         color = "Sample") +
    plot_options() +
    plot_theme() +
    theme(aspect.ratio = wide_ar)

  counts_cdf
  save_plot(counts_cdf)


}
