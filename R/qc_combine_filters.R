#' This function combines the QC filters
#' @description this function, given the zerocount filter and the low plasmid CPM filter, finds the union of the filters to produce a combined filter. It also computes and reports several other set theory descriptions of the filters
#' @param filter_zerocount A vector of TRUE/FALSE values from the zerocount filter where TRUE means that the construct was flagged by the filter
#' @param filter_plasmid A vector of TRUE/FALSE values from the low plasmid CPM filter where TRUE means that the construct was flagged by the filter
#' @param use_combined default is TRUE; whether to use a combination or the union of the two filters to produce the final QC filter
#' @importFrom magrittr %>% %<>%
#' @import dplyr
#' @return a named list with the `combined_filter` where a column `keep_pgRNA` is the final QC filter TRUE meaning that a construct passed the filters and should remain in the analysis while FALSE means the construct was flagged by filter(s) and should be removed prior to analysis
#'                      with `num_filtered_report` which reports how many and percent of constructs which were flagged by the union of the filters
#'                      with `num_which_filter_report` which reports how many and percent of constructs were flagged by only the zerocount filter, only the low plasmid cpm filter, by the intersection of the filters, and by no filters
#'

qc_combine_filters <- function(filter_zerocount,
                               filter_plasmid,
                               use_combined = TRUE) {
  combined_filter <- cbind(data.frame(filter_zerocount), data.frame(filter_plasmid)) %>%
    mutate(combined = filter_zerocount == TRUE | filter_plasmid == TRUE)

  combined_filter_df <- data.frame(
    "Filtered" = c("Yes", "No"), # Yes filtered means remove pgRNA
    n = c(sum(combined_filter$combined), sum(!combined_filter$combined))
  ) %>%
    mutate(percent = round(((n / sum(n)) * 100), 2))

  which_filter_df <- data.frame(
    "Which Filter" = c("Low plasmid cpm only", "Zero count for any sample", "Both", "Neither"),
    n = c(
      sum(combined_filter$filter_zerocount == FALSE & combined_filter$filter_plasmid == TRUE),
      sum(combined_filter$filter_zerocount == TRUE & combined_filter$filter_plasmid == FALSE),
      sum(combined_filter$filter_zerocount == TRUE & combined_filter$filter_plasmid == TRUE),
      sum(combined_filter$filter_zerocount == FALSE & combined_filter$filter_plasmid == FALSE)
    )
  ) %>%
    mutate(percent = round(((n / sum(n)) * 100), 2))

  if (use_combined) {
    combined_filter %<>% select(combined) %>%
      mutate(combined = as.logical(abs(combined - 1))) %>%
      `colnames<-`(c("keep_pgRNA")) # invert the bits so boolean reports which ones you keep, not which ones you remove.
  } else {
    combined_filter$keep_pgRNA <- NA
  }

  return(list(
    combined_filter = combined_filter$keep_pgRNA,
    num_filtered_report = combined_filter_df,
    num_which_filter_report = which_filter_df
  ))
}
