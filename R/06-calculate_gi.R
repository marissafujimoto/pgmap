#' Calculate Genetic Interaction scores
#' @description Create results table that has CRISPR scores, Wilcoxon rank-sum test and t tests.
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param test options include 'wilcoxon' and 't-test'. By default, both will be run.
#' @param overwrite default is FALSE; whether to overwrite the QC Report file
#' @param output_file default is `GI_Results`; name of the output GI results file
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
#'     replicates = "rep") %>%
#'   calc_crispr() %>%
#'   calc_gi()
#'
#' }
calc_gi <- function(.data = NULL,
                    gimap_dataset) {

  # Code adapted from
  # https://github.com/FredHutch/GI_mapping/blob/main/workflow/scripts/04-calculate_GI_scores.Rmd

  if (!is.null(.data)) gimap_dataset <- .data

  if (!("gimap_dataset" %in% class(gimap_dataset))) stop("This function only works with gimap_dataset objects which can be made with the setup_data() function.")


  ## calculate expected double-targeting GI score by summing the two mean single-targeting
  ## CRISPR scores for that paralog pair
  gi_calc_df <- gimap_dataset$crispr_score %>%
    dplyr::mutate(expected_crispr = dplyr::case_when(
                    target_type == "gene_gene" ~ mean_single_target_crispr_1 + mean_single_target_crispr_2,
                    target_type == "gene_ctrl" ~ mean_single_target_crispr_1 + mean_double_control_crispr_2,
                    target_type == "ctrl_gene" ~ mean_double_control_crispr_1 + mean_single_target_crispr_2)
                    ) %>%
    tidyr::pivot_longer(cols = c(pgRNA_target_1, pgRNA_target_2, -expected_crispr),
                        names_to = "position",
                        values_to = "pgRNA_target") %>%
    tidyr::pivot_longer(cols = c(mean_single_target_crispr_1, mean_single_target_crispr_2, -expected_crispr),
                        names_to = "target",
                        values_to = "target_crispr")

  results <- gi_calc_df %>%
    dplyr::group_by(rep, pgRNA_target) %>%
    dplyr::summarize(
      mean_expected_crispr = mean(expected_crispr, na.rm = TRUE),
      mean_observed_crispr = mean(target_crispr, na.rm = TRUE)) %>%
      group_modify(~ broom::tidy(lm(mean_observed_crispr ~ mean_expected_crispr, data = .x)))

  gimap_dataset$results <- results

  stats <- results %>%
      dplyr::ungroup() %>%
      dplyr::select(term, estimate, rep) %>%
      pivot_wider(names_from = term,
                  values_from = estimate) %>%
      rename(intercept = "(Intercept)", slope = mean_expected_crispr)

  gi_calc_adj <- gi_calc_df %>%
    left_join(stats, by = "rep") %>%
    mutate(gi_score = target_crispr - (intercept + slope * expected_crispr))

  gimap_dataset$gi_scores <- gi_calc_adj

  return(gimap_dataset)
}


#' Do tests for each sample
#' @description Create results table that has CRISPR scores, Wilcoxon rank-sum test and t tests.
#' @param .data Data can be piped in with %>% or |> from function to function. But the data must still be a gimap_dataset
#' @param gimap_dataset A special dataset structure that is setup using the `setup_data()` function.
#' @param test options include 'wilcoxon' and 't-test'. By default, both will be run.
#' @param overwrite default is FALSE; whether to overwrite the QC Report file
#' @param output_file default is `GI_Results`; name of the output GI results file
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
#'     replicates = "rep") %>%
#'   calc_crispr()
#'
#'
#'
#' }


gimap_stats <- function() {
  ## get a vector of GI scores for all single-targeting ("control") pgRNAs for each rep
  ## get double-targeting pgRNAs for this rep, do a t-test to compare the double-

  ## targeting GI scores for each paralog pair to the control vector

  ## adjust for multiple testing using the Benjamini-Hochberg method

  d.double_GI_scores <- gi_calc_adj %>%
    filter(rep == i) %>%
    group_by(paralog_pair) %>%
    mutate(p_val = t.test(x = single_GI_scores,
                          y = GI_score,
                          paired = FALSE)$p.value)

  ## adjust for multiple testing using the Benjamini-Hochberg method
  d.p_val <- d.double_GI_scores %>%
    dplyr::select(paralog_pair, p_val) %>%
    arrange(p_val) %>%
    distinct(p_val, .keep_all = TRUE)

  p_vals <- d.p_val %>%
    pull(p_val)

  fdr_vals <- p.adjust(p_vals, method = "BH")

  fdr_df <- tibble("fdr" = fdr_vals) %>%
    bind_cols(d.p_val) %>%
    dplyr::select(-p_val)

  ## add FDR values back into the double-targeting DF
  d.double_GI_scores <- left_join(d.double_GI_scores, fdr_df, by = "paralog_pair")

}
