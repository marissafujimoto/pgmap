#' Normalize Log fold changes
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
#'   gimap_normalize(
#'     timepoints = "day",
#'     replicates = "rep")
#'
#' # To see results
#' gimap_dataset$normalized_log_fc
#' }
gimap_normalize <- function(.data = NULL,
                     gimap_dataset,
                     timepoints = NULL,
                     replicates = NULL) {

  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

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

  message("Normalizing Log Fold Change")

  # Based on log fold change calculations and other handling will go based on the code in:
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/03-filter_and_calculate_LFC.Rmd

  # Doing some reshaping to get one handy data frame
  lfc_df <- gimap_dataset$transformed_data$log2_cpm %>%
    as.data.frame() %>%
    dplyr::mutate(pg_ids = gimap_dataset$metadata$pg_ids$id) %>%
    tidyr::pivot_longer(-pg_ids) %>%
    dplyr::left_join(gimap_dataset$metadata$sample_metadata, by = c("name" = "col_names")) %>%
    dplyr::select(-name) %>%
    # Take an average for
    tidyr::pivot_wider(values_from = value,
                       names_from = c(timepoints, replicates))

  late_df <- dplyr::select(lfc_df, dplyr::starts_with("late"))
  early_df <- apply(dplyr::select(lfc_df, dplyr::starts_with("early")), 1, mean)
  plasmid_df <- apply(dplyr::select(lfc_df, starts_with("plasmid")), 1, mean)


  # Do the actual calculations
  #late_vs_early <-  lfc_df %>%
  #  dplyr::mutate_at(dplyr::vars(dplyr::starts_with("late")), ~.x - early_df ) %>%
  #  dplyr::select(pg_ids, dplyr::starts_with("late"))

  late_vs_plasmid_df <-  lfc_df %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("late")), ~.x - plasmid_df) %>%
    dplyr::select(pg_ids, dplyr::starts_with("late"))  %>%
    dplyr::left_join(gimap_dataset$annotation, by = c("pg_ids" = "pgRNA_id"))


  ########################### Perform adjustments #############################

  ### Calculate medians
  neg_control_median_df <- late_vs_plasmid_df %>%
    dplyr::filter(norm_ctrl_flag == "negative_control") %>%
    dplyr::select(pg_ids, dplyr::starts_with("late"))

  neg_control_median <- median(apply(neg_control_median_df[, -1], 1, median))

  # First and second adjustments to LFC
  lfc_df <- lfc_df %>%
    tidyr::pivot_longer(dplyr::starts_with("late"),
                        names_to = "rep",
                        values_to = "lfc_plasmid_vs_late")  %>%
    dplyr::left_join(gimap_dataset$annotation, by = c("pg_ids" = "pgRNA_id")) %>%
    dplyr::mutate(
      # Take LFC, then subtract median of negative controls.
      # This will result in the median of the nontargeting being set to 0.
      lfc_adj = lfc_plasmid_vs_late - neg_control_median,
      # Then, divide by the median of negative controls (double non-targeting) minus
      # median of positive controls (targeting 1 essential gene).
      # This will effectively set the median of the positive controls (essential genes) to -1.
      lfc_adj = lfc_adj / (median(lfc_adj[norm_ctrl_flag == "negative_control"]) - median(lfc_adj[norm_ctrl_flag == "positive_control"]))
    )

  # Save this at the construct level
  gimap_dataset$normalized_log_fc <- lfc_df

  return(gimap_dataset)
}
