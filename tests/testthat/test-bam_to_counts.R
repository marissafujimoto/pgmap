test_that("bam to counts conversion", {

   bam_dir <- file.path(example_data_folder(), "bam") 
   
   sample_names <- c("sample1", "sample2", "sample3")
   
   counts <- calc_counts(
     bam_dir = bam_dir, 
     sample_names = sample_names
  )
   
   # Expect the names with the samples 
   testthat::expect_named(counts$counts, c("id", "sample1", "sample2", "sample3"))
   testthat::expect_named(counts$stats, c("sample", "n_correctly_paired", "n_total"))
   
   # If the number of ids change we wanna know
   testthat::expect_length(counts$counts$id, 32508) 
   
   
   # This is just a test to see if the numbers changed
   subset_data <- as.matrix(counts$counts[1, -1])
   dimnames(subset_data) <- NULL
   test <- matrix(c(5,6,4), nrow = 1, ncol = 3)
   testthat::expect_equal(subset_data,test)
})
