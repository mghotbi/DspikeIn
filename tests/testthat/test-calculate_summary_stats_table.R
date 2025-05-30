library(testthat)
library(DspikeIn)
library(flextable)
library(dplyr)
library(utils)
library(microbiome)

test_that("calculate_summary_stats_table() returns a flextable and saves files", {
  data("physeq_16SOTU", package = "DspikeIn")
  absolute_count <- microbiome::meta(physeq_16SOTU)

  expect_true(nrow(absolute_count) > 0, info = "Dataset should have rows.")
  expect_true(ncol(absolute_count) > 0, info = "Dataset should have columns.")

  temp_docx <- tempfile(fileext = ".docx")
  temp_csv <- sub(".docx", ".csv", temp_docx)

  summary_table <- calculate_summary_stats_table(absolute_count, output_path = temp_docx)

  expect_s3_class(summary_table, "flextable")
  expect_true(file.exists(temp_docx), info = "Word document should be created.")
  expect_true(file.exists(temp_csv), info = "CSV file should be created.")

  summary_data <- read.csv(temp_csv)
  expect_true(nrow(summary_data) > 0, info = "CSV file should have rows.")
  expect_true(ncol(summary_data) > 0, info = "CSV file should have columns.")

  # Cleanup
  unlink(c(temp_docx, temp_csv))
})
