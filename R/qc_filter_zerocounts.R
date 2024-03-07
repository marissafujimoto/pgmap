#' Create a filter for pgRNAs which have a raw count of 0 for any sample/time point
#' @description This function flags and reports which and how many pgRNAs have a raw count of 0 for any sample/time point
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @return a named list with the filter `filter` specifying which pgRNA have a count zero for at least one sample/time point and a report df `reportdf` for the number and percent of pgRNA which have a count zero for at least one sample/time point
#' @examples \dontrun{
#'
#' }
#'
qc_filter_zerocounts <- function(gimap_dataset) {

  counts_filter <- unlist(lapply(1:nrow(gimap_dataset$raw_counts), function(x) 0 %in% gimap_dataset$raw_counts[x, ]))

  zerocount_df <- data.frame("RawCount0" = c(FALSE, TRUE), n = c(sum(!counts_filter), sum(counts_filter))) %>%
    mutate(percent = round(((n / sum(n)) * 100), 2))

  return(list(filter = counts_filter, reportdf = zerocount_df))
}
