#'
#'
#'
#'
#'
#' @importFrom magrittr %>% %<>%
#'

qc_combine_filters <- function(filter_zerocount, filter_plasmid, use_combined = TRUE){
  combined_filter <- cbind(data.frame(filter_zerocount), data.frame(filter_plasmid)) %>%
    mutate(combined = filter_zerocount == TRUE | filter_plasmid == TRUE)
  
  combined_filter_df <- data.frame("Filtered" = c("Yes", "No"), #Yes filtered means remove pgRNA
                                   n = c(sum(combined_filter$combined), sum(!combined_filter$combined))) %>%
    mutate(percent = round(((n/sum(n))*100),2))
  
  which_filter_df <- data.frame("Which Filter" = c("Low plasmid cpm only", "Zero count for any sample", "Both", "Neither"),
                                 n = c(sum(combined_filter$filter_zerocount == FALSE & combined_filter$filter_plasmid == TRUE),
                                              sum(combined_filter$filter_zerocount == TRUE & combined_filter$filter_plasmid == FALSE),
                                              sum(combined_filter$filter_zerocount == TRUE & combined_filter$filter_plasmid == TRUE),
                                              sum(combined_filter$filter_zerocount == FALSE & combined_filter$filter_plasmid == FALSE))) %>%
    mutate(percent = round(((n/sum(n)) *100),2))
  
  if (use_combined){
    combined_filter %<>% select(combined) %>% mutate(combined = as.logical(abs(combined -1))) %>% `colnames<-`(c("keep_pgRNA")) #invert the bits so boolean reports which ones you keep, not which ones you remove.
  } else { combined_filter$keep_pgRNA <- NA }
  
  return(list(combined_filter = combined_filter, 
              num_filtered_report = combined_filter_df,
              num_which_filter_report = which_filter_df))
  
}