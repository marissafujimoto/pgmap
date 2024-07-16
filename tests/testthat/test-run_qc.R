test_that("HTML file is created and content is correct", {
  gimap_dataset <- get_example_data("gimap")
  html_file <- run_qc(gimap_dataset)

  expect_true(file.exists(html_file))
})
