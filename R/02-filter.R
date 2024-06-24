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
#' 

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
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom purrr reduce
#' @return a named list with the filter `filter` specifying which pgRNA have a count zero for at least one sample/time point and a report df `reportdf` for the number and percent of pgRNA which have a count zero for at least one sample/time point 
#' @examples \dontrun{
#' 
#' }
#'

qc_filter_zerocounts <- function(gimap_dataset){
  
  counts_filter <- data.frame(gimap_dataset$raw_counts) %>% map(~.x %in% c(0)) %>% reduce(`|`)
  
  zerocount_df <- data.frame("RawCount0" = c(FALSE, TRUE), n = c(sum(!counts_filter), sum(counts_filter))) %>%
    mutate(percent = round(((n/sum(n))*100),2))
  
  return(list(filter = counts_filter, reportdf = zerocount_df))
  
}