#' Returns example data for package
#' @description This function loads and returns example data for the packagae. Which dataset is returned must be specified
#' @param which_data options are "count" or "meta"; specifies which example dataset should be returned
#' @export
#' @examples \dontrun{
#'
#' pg_data <- get_example_data()
#' }
get_example_data <- function(which_data) {
  if (which_data == "count") {
    file <- list.files(
      pattern = "PP_pgPEN_HeLa_counts.txt",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_tsv(file))
  } else if (which_data == "meta") {
    file <- list.files(
      pattern = "pgRNA_ID_pgPEN_library_comp.csv",
      recursive = TRUE,
      system.file("extdata", package = "gimap"),
      full.names = TRUE
    )
    return(readr::read_csv(file, skip = 1))
  } else {
    stop("Specification for `which_data` not understood; Need to use 'count' or 'meta'")
  }
}

#' @import ggplot2

## ggplot themes
## see: https://www.rdocumentation.org/packages/ggplot2/versions/2.1.0/topics/theme_update
## and https://stackoverflow.com/questions/23173915/can-ggplot-theme-formatting-be-saved-as-an-object
plot_theme <- function() {
  theme(
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(color = "black")
  )
}

#' @import ggplot2
plot_options <- function() {
  list(theme_bw(base_size = 12))
}

#' @import kableExtra
print_kbl <- function(tbl) {
  kbl(tbl) %>%
    kable_styling(
      full_width = FALSE,
      position = "left",
      bootstrap_options = c("striped", "hover", "responsive")
    )
}

save_tbl <- function(tbl, out_dir = NULL, params = NULL) {
  tbl_str <- deparse(substitute(tbl))
  tbl_name <- str_split(tbl_str, pattern = "\\.")[[1]][2]
  write_tsv(tbl, file.path(out_dir, "tables", "tsv", paste0(params$cell_line, "_", tbl_name, ".txt")))
  write_rds(tbl, file.path(out_dir, "tables", "rds", paste0("d.", params$cell_line, "_", tbl_name)))
}

#' @import ggplot2
save_plot <- function(plt, out_dir = NULL, params = list(cell_line = NULL)) {
  plt_str <- deparse(substitute(plt))
  if (!dir.exists(file.path(out_dir, "plots", "pdf"))) {
    dir.create(file.path(out_dir, "plots", "pdf"), recursive = TRUE)
  }
  ggsave(
    plot = plt,
    filename = file.path(out_dir, "plots", "pdf", paste0(params$cell_line, "_", plt_str, ".pdf"))
  )
  if (!dir.exists(file.path(out_dir, "plots", "png"))) {
    dir.create(file.path(out_dir, "plots", "png"), recursive = TRUE)
  }
  ggsave(
    plot = plt,
    filename = file.path(out_dir, "plots", "png", paste0(params$cell_line, "_", plt_str, ".png"))
  )
}

make_out_dir <- function(out_dir) {
  if (dir.exists(out_dir)) {
    ## print a message with the output directory location
    print(paste("Output directory already exists at:", out_dir, sep = " "))
  } else {
    ## make output dirs
    dir.create(file.path(out_dir, "tables", "rds"), recursive = TRUE)
    dir.create(file.path(out_dir, "tables", "tsv"), recursive = TRUE)
    dir.create(file.path(out_dir, "plots", "png"), recursive = TRUE)
    dir.create(file.path(out_dir, "plots", "pdf"), recursive = TRUE)
    print(paste("Output directory created at:", out_dir, sep = " "))
  }
}
