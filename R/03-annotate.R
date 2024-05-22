#' Annotate gimap data
#' @description In this function, a `gimap_dataset` is annotated as far as which genes should be used as controls.
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param cell_line which cell line are you using? Default is "HELA"
#' @param cn_annotate TRUE or FALSE you'd also like to have Copy number annotation from DepMap. These data are optional
#' @param annotation_file If no file is given, will attempt to use the design file from https://media.addgene.org/cms/filer_public/a9/9a/a99a9328-324b-42ff-8ccc-30c544b899e4/pgrna_library.xlsx
#' @param control_genes A vector of gene symbols (e.g. AAMP) that should be labeled as control genes. These will be used for log fold change calculations. If no list is given then DepMap Public 23Q4 Achilles_common_essentials.csv is used https://depmap.org/portal/download/all/
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
#' gimap_dataset$annotation
#' }
gimap_annotate <- function(.data = NULL,
                           gimap_dataset,
                           cell_line = "HELA",
                           control_genes = NULL,
                           cn_annotate = TRUE,
                           annotation_file = NULL) {

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  # Get the annotation data based on the pg construct design
  if (!is.null(annotation_file)) {
    if (!file.exists(annotation_file)) stop("The annotation_file specified cannot be found. Please double check the file path")
    annotation_df <- read_table(annotation_file)
  } else {
    annotation_df <- get_example_data("annotation")
  }

  ############################ CONTROL GENE ANNOTATION #########################
  # If control genes aren't provided then we get some from DepMap
  if (!is.null(control_genes)) {
    if (!file.exists(control_genes)) stop("The annotation_file specified cannot be found. Please double check the file path")
    control_genes <- read_table(control_genes)[, 1]
  } else {
    # This file is from https://depmap.org/portal/download/all/ and from DepMap Public 19Q3 All Files
    # Essential gene labeling is from inst/extdata/Achilles_common_essentials.csv
    control_genes <- readr::read_tsv("https://figshare.com/ndownloader/files/40448429", show_col_types = FALSE)
    control_genes <- control_genes %>%
      tidyr::separate(col = Gene, into = c("gene_symbol", "entrez_id"), remove = FALSE, extra = "drop")
  }

  ############################ Get TPM data ####################################
  # This is not optional because its used to flag things
  ## get TPM and CN information (w/ option for user to upload their own info)
  depmap_metadata <- readr::read_csv("https://figshare.com/ndownloader/files/35020903", show_col_types = FALSE)

  my_depmap_id <- depmap_metadata %>%
    dplyr::filter(stripped_cell_line_name == cell_line) %>%
    dplyr::pull(DepMap_ID)

  tpm_file <- file.path(system.file("extdata", package = "gimap"), "CCLE_expression.csv")

  if (!file.exists(tpm_file)) tpm_setup()

  depmap_tpm <- readr::read_csv(tpm_file,
    show_col_types = FALSE,
    col_select = c("genes", my_depmap_id)
  ) %>%
    dplyr::rename(log2_tpm = my_depmap_id) %>%
    dplyr::mutate(expressed_flag = dplyr::case_when(
      log2_tpm < 1 ~ FALSE,
      log2_tpm >= 1 ~ TRUE,
      is.na(log2_tpm) ~ NA
    ))

  ############################ COPY NUMBER ANNOTATION ##########################
  if (cn_annotate) {

    cn_file <- file.path(system.file("extdata", package = "gimap"), "CCLE_gene_cn.csv")
    if (!file.exists(cn_file)) cn_setup()

    # Read in the CN data
    depmap_cn <- readr::read_csv(cn_file,
      show_col_types = FALSE,
      col_select = c("genes", my_depmap_id)
    ) %>%
      dplyr::rename(log2_cn = my_depmap_id)

    annotation_df <- annotation_df %>%
      dplyr::left_join(depmap_cn, by = c("gene1_symbol" = "genes")) %>%
      dplyr::left_join(depmap_cn, by = c("gene2_symbol" = "genes"), suffix = c("_gene1", "_gene2"))
  }

  ############################ ANNOTATION COMBINING ############################
  # This set up is more or less the same as the original
  # https://github.com/FredHutch/GI_mapping/blob/41ac7d5ed7025252343e2c823fba22f8a363e25c/workflow/scripts/02-get_pgRNA_annotations.Rmd#L435
  annotation_df <- annotation_df %>%
    dplyr::left_join(depmap_tpm, by = c("gene1_symbol" = "genes")) %>%
    dplyr::rename(gene1_essential_flag = expressed_flag) %>%
    dplyr::left_join(depmap_tpm, by = c("gene2_symbol" = "genes"), suffix = c("_gene1", "_gene2")) %>%
    dplyr::rename(gene2_essential_flag = expressed_flag) %>%
    dplyr::mutate(norm_ctrl_flag = dplyr::case_when(
      target_type == "gene_gene" ~ "double_targeting",
      target_type == "gene_ctrl" & gene1_essential_flag == TRUE ~ "positive_control",
      target_type == "ctrl_gene" & gene2_essential_flag == TRUE ~ "positive_control",
      target_type == "gene_ctrl" & gene1_essential_flag != TRUE ~ "single_targeting",
      target_type == "ctrl_gene" & gene2_essential_flag != TRUE ~ "single_targeting",
      target_type == "ctrl_ctrl" ~ "negative_control"
    )) %>%
    dplyr::mutate(norm_ctrl_flag = factor(norm_ctrl_flag, levels = c(
      "negative_control",
      "positive_control",
      "single_targeting",
      "double_targeting"
    )))

  ################################ STORE IT ####################################
  gimap_dataset$annotation <- annotation_df

  return(gimap_dataset)
}
