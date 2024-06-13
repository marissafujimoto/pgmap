#' A function to run filtering
#' @description This is a function here's where we describe what it does
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param filter_type Can be one of the following: `zero_count_only`, `low_plasmid_cpm_only` or `rep_variation`, `zero_in_last_time_point` or a vector that includes multiple of these filters. 
#' You should decide on the appropriate filter based on the results of your QC report.
#' @returns a filtered version of the gimap_dataset returned in the $filtered_data section
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#'
#' # Highly recommended but not required
#' run_qc(gimap_dataset)
#'
#'
#' gimap_dataset <- gimap_filter(gimap_dataset)
#'
#' # To see filtered data
#' gimap_dataset$filtered_data
#'
#' }
gimap_filter <- function(.data = NULL,
                         gimap_dataset,
                         filter_type = "both") {

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  gimap_dataset$filtered <- NULL #TODO: Filtered version of the data can be stored here

  return(gimap_dataset)
}

# Possible filters

#' Create a filter for pgRNAs which have a raw count of 0 for any sample/time point
#' @description This function flags and reports which and how many pgRNAs have a raw count of 0 for any sample/time point
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the raw count data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @return a named list with the filter `filter` specifying which pgRNA have a count zero for at least one sample/time point and a report df `reportdf` for the number and percent of pgRNA which have a count zero for at least one sample/time point 
#' @examples \dontrun{
#'   qc_filter_zerocounts(gimap_dataset)
#' }
#'

qc_filter_zerocounts <- function(gimap_dataset){
  
  counts_filter <- unlist(lapply(1:nrow(gimap_dataset$raw_counts), function(x) 0 %in% gimap_dataset$raw_counts[x,]))
  
  zerocount_df <- data.frame("RawCount0" = c(FALSE, TRUE), n = c(sum(!counts_filter), sum(counts_filter))) %>%
    mutate(percent = round(((n/sum(n))*100),2))
  
  return(list(filter = counts_filter, reportdf = zerocount_df))
  
}

#' Create a filter for pgRNAs which have a low log2 CPM value for the plasmid/Day 0 sample/time point
#' @description This function flags and reports which and how many pgRNAs have low log2 CPM values for the plasmid/Day 0 sample/time point
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the log2 CPM transformed data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @return a named list with the filter `plasmid_filter` specifying which pgRNAs have low plasmid log2 CPM and a report df `plasmid_filter_report` for the number and percent of pgRNA which have a low plasmid log2 CPM
#' @examples \dontrun{
#'   qc_filter_plasmid(gimap_dataset)
#'   
#'   qc_filter_plasmid(gimap_dataset, cutoff=2)
#' }
#'

qc_filter_plasmid <- function(gimap_dataset, cutoff = NULL){
  plasmid_data <- data.frame(gimap_dataset$transformed_data$log2_cpm[, 1]) %>% `colnames<-`(c("log2_cpm"))
  
  if (is.null(cutoff)) {
    # if cutoff is null, use lower outlier 
    quantile_info <- quantile(plasmid_data$log2_cpm)
    
    cutoff <- quantile_info["25%"] - (1.5 * (quantile_info["75%"] - quantile_info["25%"])) #later step make a function for this in utils since it's used more than once
  }
  
  plasmid_cpm_filter <- unlist(lapply(1:nrow(plasmid_data), function(x) plasmid_data$log2_cpm[x] < cutoff))
  
  plasmid_filter_df <- data.frame("Plasmid_log2cpmBelowCutoff" = c(FALSE, TRUE), n = c(sum(!plasmid_cpm_filter), sum(plasmid_cpm_filter))) %>%
    mutate(percent = round(((n / sum(n)) * 100), 2)) #later step make a function for this in utils since it's used more than once
  
  return(list(
    plasmid_filter = plasmid_cpm_filter,
    plasmid_filter_report = plasmid_filter_df
  ))
}