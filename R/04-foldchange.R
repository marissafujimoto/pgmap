#' Calculate log fold change for a
#' @description This calculates the log fold change for a gimap dataset based on the annotation and metadata provided.
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param timepoints Specifies the column name of the metadata set up in `$metadata$sample_metadata` that has a factor that represents the timepoints. The column used for timepoints must be numeric or at least ordinal.
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
#'   gimap_annotate() %>%
#'   calc_lfc()
#'
#' # To see results
#' gimap_dataset$log_fc
#' }
calc_lfc <- function(.data = NULL,
                     gimap_dataset,
                     timepoints = NULL,
                     replicates = NULL) {

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")

  if (!is.null(replicates)) {
    if (!(replicates %in% colnames(gimap_dataset$metadata$sample_metadata))) {
      stop("The column name specified for 'replicates' does not exist in gimap_dataset$metadata$sample_metadata")
    }
    # Rename and recode the timepoints variable
    gimap_dataset$metadata$sample_metadata <- gimap_dataset$metadata$sample_metadata %>%
      dplyr::rename(replicates = all_of(replicates))
  }

  if (!is.null(timepoints)) {
    if (!(timepoints %in% colnames(gimap_dataset$metadata$sample_metadata))) {
      stop("The column name specified for 'timepoints' does not exist in gimap_dataset$metadata$sample_metadata")
    }

    # Rename and recode the timepoints variable
    gimap_dataset$metadata$sample_metadata <- gimap_dataset$metadata$sample_metadata %>%
      dplyr::rename(timepoints = all_of(timepoints)) %>%
      dplyr::mutate(timepoints = dplyr::case_when(
        timepoints == min(timepoints) ~ "plasmid",
        timepoints == max(timepoints) ~ "late",
        TRUE ~ "early"))
  }

  if (is.null(gimap_dataset$annotation)) {
    stop(
      "No annotations are stored in this gimap_dataset, annotations are needed so we know what genes should be used as controls.",
      "Please run gimap_annotate() function on your gimap_dataset and then retry this function."
    )
  }

  message("Calculating Log Fold Change")

  # Based on log fold change calculations and other handling will go based on the code in:
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

  # Doing some reshaping to get one handy data frame
  lfc_df <- gimap_dataset$transformed_data$log2_cpm %>%
    as.data.frame() %>%
    dplyr::mutate(pg_ids = gimap_dataset$metadata$pg_ids$id) %>%
    tidyr::pivot_longer(-pg_ids) %>%
    dplyr::left_join(gimap_dataset$metadata$sample_metadata, by = c("name" = "col_names")) %>%
    dplyr::group_by(timepoints, pg_ids) %>%
    dplyr::summarize(timepoint_avg = mean(value)) %>%
    tidyr::pivot_wider(values_from = timepoint_avg,
                       names_from = timepoints)

  # Do the actual calculations
  lfc_df$lfc_plasmid_vs_late <- lfc_df$late - lfc_df$plasmid
  lfc_df$lfc_early_vs_late <- lfc_df$late - lfc_df$early

  ########################### Perform adjustments #############################

  lfc_df <- lfc_df %>%
    dplyr::left_join(gimap_dataset$annotation, by = c("pg_ids" = "pgRNA_id"))

  ### Calculate medians
  neg_control_median <- median(lfc_df$lfc_plasmid_vs_late[lfc_df$norm_ctrl_flag == "negative_control"])

  # First and second adjustments to LFC
  lfc_df <- lfc_df %>%
    dplyr::mutate(
      # Take LFC, then subtract median of negative controls.
      # This will result in the median of the nontargeting being set to 0.
      lfc_adj = lfc_plasmid_vs_late - neg_control_median,
      # Then, divide by the median of negative controls (double non-targeting) minus
      # median of positive controls (targeting 1 essential gene).
      # This will effectively set the median of the positive controls (essential genes) to -1.
      lfc_adj = lfc_adj / (median(lfc_adj[norm_ctrl_flag == "negative_control"]) - median(lfc_adj[norm_ctrl_flag == "positive_control"]))
    )

  # Calculate medians based on single, double targeting as well as if they are unexpressed control genes
  medians_df <- lfc_df %>%
    dplyr::group_by(target_type, unexpressed_ctrl_flag) %>%
    dplyr::summarize(median = median(lfc_adj)) %>%
    dplyr::filter(unexpressed_ctrl_flag)

  # Third adjustment
  lfc_df <- lfc_df %>%
    dplyr::left_join(medians_df, by = "target_type") %>%
    dplyr::mutate(
      # Since the pgPEN library uses non-targeting controls, we adjusted for the
      # fact that single-targeting pgRNAs generate only two double-strand breaks
      # (1 per allele), whereas the double-targeting pgRNAs generate four DSBs.
      # To do this, we set the median (adjusted) LFC for unexpressed genes of each group to zero.
      crispr_score = dplyr::case_when(
        target_type == "single_targeting" ~ lfc_adj - median,
        target_type == "double_targeting" ~ lfc_adj - median,
        TRUE ~ lfc_adj
      ),
      n_genes_expressed = dplyr::case_when(
        gene1_expressed_flag == FALSE & gene2_expressed_flag == FALSE ~ "0",
        gene1_expressed_flag == TRUE & gene2_expressed_flag == FALSE ~ "1",
        gene1_expressed_flag == FALSE & gene2_expressed_flag == TRUE ~ "1",
        gene1_expressed_flag == TRUE & gene2_expressed_flag == TRUE ~ "2")
    ) %>%
    dplyr::select(pgRNA_target,
                  pg_ids,
                  early,
                  late,
                  plasmid,
                  lfc_plasmid_vs_late,
                  lfc_early_vs_late,
                  lfc_adj,
                  crispr_score)

  # Save this at the construct level
  gimap_dataset$log_fc <- lfc_df

  # Summarize to target level and save that
  crispr_df <- lfc_df %>%
    dplyr::group_by(pgRNA_target) %>%
    dplyr::summarize(target_mean_cs = mean(crispr_score),
           target_median_cs = median(crispr_score)) %>%
    # TODO: This is suspicious step its going to drop a lot of data
    dplyr::select(-pg_ids)

  # Save at the target level
  gimap_dataset$crispr <- crispr_df

  return(gimap_dataset)
}
