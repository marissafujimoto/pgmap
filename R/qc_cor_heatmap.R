#' Create a correlation heatmap for the pgRNA cpm's
#' @description This function uses the cor function to find correlations between the sample CPM's and then plots a heatmap of these
#' @param gimap_dataset The special gimap_dataset from the `setup_data` function which contains the transformed data
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @return sample_cor_heatmap a pheatmap
#' @examples \dontrun{
#' 
#' }
#'

qc_cor_heatmap <- function(gimap_dataset){
  
  cpm_cor <- gimap_dataset$transformed_data$cpm %>%
    cor() %>%
    round(2) %>%
    data.frame()
  
  sample_cor_heatmap <-
    pheatmap(cpm_cor,
             border_color = "white",
             cellwidth = 20, cellheight = 20,
             treeheight_row = 20, treeheight_col = 20)
  
  return (sample_cor_heatmap)
  
}