#' Obtain the counts for a group of samples
#' @description This function takes a forward and reverse bam files for a group of samples and returns the counts and stats.
#' @param bam_dir a file path to a directory where the bam files associated with the sample names will be stored.
#' @param sample_names a character vector that indicates the name of the samples as they are listed in the file names.
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @importFrom furrr future_pmap future_map
#' @importFrom purrr reduce
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @export
#' @examples \dontrun{
#'
#' # This will be the directory that contains all the file path
#' bam_dir <- file.path(example_data_folder(), "bam")
#'
#' # These MUST be names that are listed in the bam file name itself.
#' # There needs to be exactly 2 bam files per sample
#' sample_names <- c("sample1", "sample2", "sample3")
#'
#' counts_df <- calc_counts(
#'   bam_dir = bam_dir,
#'   sample_names = sample_names,
#'   time = TRUE
#' )
#'
#' # We can write these to CSV files like this
#' readr::write_csv(counts_df$counts, "counts_pgmap.tsv")
#' readr::write_csv(counts_df$stats, "stats_pgmap.tsv")
#' }
calc_counts <- function(bam_dir, sample_names, time = TRUE) {
  if (time) timing <- Sys.time()

  # Grab the file names
  sample_df <- grab_paired_files(
    dir = bam_dir,
    sample_names = sample_names
  )

  # Now use the custom function to get the counts for each sample from the pair of bam files
  counts <- furrr::future_pmap(sample_df, function(sample_name, forward_file, reverse_file) {
    sample_count(
      bam_1 = forward_file,
      bam_2 = reverse_file,
      sample_name = sample_name
    )
  })
  # Bring along sample names
  names(counts) <- sample_names

  # Then condense all samples to one stats data frame
  stats_df <- furrr::future_map(counts, "stats") %>%
    dplyr::bind_rows(.id = "sample")

  # Then condense all samples to one counts data frame
  counts_df <- furrr::future_map(counts, "counts", right_join) %>%
    purrr::reduce(dplyr::full_join)

  if (time) message(paste0("Done(", signif(as.numeric(difftime(Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))
  # Return as a list of two data frames
  return(list(counts = counts_df, stats = stats_df))
}

#' Obtain the counts for a single sample
#' @description This function takes a forward and reverse bam file and returns the counts
#' @param bam_1 a file path to one of the bam files for the sample
#' @param bam_2 a file path to the other one of the bam files for the sample
#' @param sample_name a single character that indicates the name of the sample "sample1"
#' @param time TRUE/FALSE you want the duration this takes to run printed out
#' @import dplyr
#' @importFrom Rsamtools ScanBamParam scanBam scanBamFlag
#' @export
#' @examples \dontrun{
#'
#' bam_dir <- file.path(example_data_folder(), "bam")
#'
#' counts <- sample_count(
#'   bam_1 = file.path(bam_dir, "pgMAP_tutorial_gRNA1_trimmed_sample1_aligned.bam"),
#'   bam_2 = file.path(bam_dir, "pgMAP_tutorial_gRNA2_trimmed_sample1_aligned.bam"),
#'   sample_name = "sample1",
#'   time = TRUE
#' )
#' }
sample_count <- function(bam_1, bam_2, sample_name, time = FALSE) {
  if (time) timing <- Sys.time()
  # Print out message
  message(paste("Parsing reads for", sample_name))

  # define parameters for reading in BAM files, then read them in:
  # qname = the name of the mapped read (query template name)
  # rname = the name of the reference sequence that the read aligned to (reference sequence name)
  param <- ScanBamParam(
    ## restrict to mapped reads
    flag = scanBamFlag(isUnmappedQuery = FALSE),
    ## only read in the necessary fields
    what = c("qname", "rname")
  )

  # Read in the BAM files
  bam_1 <- scanBam(bam_1, param = param) %>%
    as.data.frame()
  bam_2 <- scanBam(bam_2, param = param) %>%
    as.data.frame()

  # Join forward and reverse bams together.
  paired_df <- dplyr::inner_join(bam_1, bam_2,
    by = "qname",
    relationship = "many-to-many",
    suffix = c("_1", "_2")
  ) %>%
    dplyr::mutate(paired = rname_1 == rname_2)

  # If a given set of reads have one or more correct pairings, then keep the
  # correct pairings and discard all incorrect pairings for those reads.
  qname2anypaired <- paired_df %>%
    dplyr::group_by(qname) %>%
    dplyr::summarize(any_paired = as.logical(any(paired))) %>%
    dplyr::collect()

  # Calculate weights
  # This contains code from the original pipeline that I don't really understand but does seem to effect final numbers
  weights_df <- paired_df %>%
    dplyr::left_join(qname2anypaired, by = "qname") %>%
    filter(paired | !any_paired) %>%
    select(-any_paired) %>%
    dplyr::group_by(qname) %>%
    dplyr::count() %>%
    dplyr::mutate(weight = 1 / n)

  # Join weights to original data
  paired_df <- paired_df %>%
    dplyr::left_join(weights_df, by = c("qname")) %>%
    dplyr::select(qname, rname = rname_1, n, weight, paired)

  # Getting the counts
  counts_df <- paired_df %>%
    dplyr::filter(paired) %>%
    dplyr::select("id" = rname, weight) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(sample_name = sum(weight)) %>%
    dplyr::ungroup()

  # rename so sample id is in the column name
  colnames(counts_df) <- c("id", sample_name)

  # Calculate stats
  n_paired <- paired_df %>%
    dplyr::filter(paired) %>%
    dplyr::summarize(sum(weight)) %>%
    dplyr::collect() %>%
    .[[1]]
  n <- paired_df %>%
    dplyr::summarize(sum(weight)) %>%
    dplyr::collect() %>%
    .[[1]]
  stats_df <- tibble("n_correctly_paired" = n_paired, "n_total" = n)

  if (time) message(paste0("Done(", signif(as.numeric(difftime(Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))
  return(list(counts = counts_df, stats = stats_df))
}
