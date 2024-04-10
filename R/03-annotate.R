#' Annotate gimap data
#' @description In this function, a `gimap_dataset` is annotated as far as which genes should be used as controls.
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param annotation_file If no file is given, will attempt to use the design file from https://media.addgene.org/cms/filer_public/a9/9a/a99a9328-324b-42ff-8ccc-30c544b899e4/pgrna_library.xlsx
#' @param control_genes A list of genes that should be labeled as control genes. These will be used for log fold change calculations. If no list is given then DepMap Public 23Q4 Achilles_common_essentials.csv is used https://depmap.org/portal/download/all/
#' @export
#' @examples \dontrun{
#'
#' gimap_dataset <- get_example_data("gimap")
#'
#' # Highly recommended but not required
#' run_qc(gimap_dataset)
#'
#' gimap_dataset <- gimap_dataset %>%
#'   gimap_filter() %>%
#'   gimap_annotate()
#'
#' # To see anotations
#' gimap_dataset$annotations
#' }
gimap_annotate <- function(.data = NULL,
                           gimap_dataset,
                           annotation_file = NULL) {

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  if (!is.null(annotation_file)) {
    if (!file.exists(annotation_file)) stop("The annotation_file specified cannot be found. Please double check the file path")
    annotation_df <- read_table(annotation_file)
  } else {
    annotation_df <- get_example_data("annotation")
  }

  if (!is.null(control_genes)) {
    if (!file.exists(control_genes)) stop("The annotation_file specified cannot be found. Please double check the file path")
    control_genes <- read_table(control_genes)[, 1]
  } else {
    # This file is from https://depmap.org/portal/download/all/ and from DepMap Public 19Q3 All Files
    # Essential gene labeling is from inst/extdata/Achilles_common_essentials.csv
    control_genes <- get_example_data("crtl_genes")
    control_genes <- stringr::word(control_genes$gene, sep = " \\(", 1)
  }

  annotation_df <-annotation_df %>%
    dplyr::mutate(
      gene_essential_flag = gene1_symbol %in% control_genes | gene2_symbol %in% control_genes
      )

  # This is the core of the code we'll need but we'll need to refactor
  annotation_df  <- annotation_df  %>%
  dplyr::mutate(norm_ctrl_flag = dplyr::case_when(
    target_type == "gene_gene" ~ "double_targeting",
    gene_essential_flag == TRUE ~ "positive_control",
    gene_essential_flag != TRUE ~ "single_targeting",
    target_type == "ctrl_ctrl" ~ "negative_control")) %>%
  dplyr::mutate(norm_ctrl_flag = factor(norm_ctrl_flag, levels = c("negative_control",
                                                            "positive_control",
                                                            "single_targeting",
                                                            "double_targeting")))

  gimap_dataset$annotation <- NULL #TODO: Final step is annotations that line up to the same order as the pg gene data should be stored here.

  return(gimap_dataset)
}
