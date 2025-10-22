library(testthat)
library(DspikeIn)
library(dplyr)
library(utils)
library(microbiome)

test_that("calculate_summary_stats_table() saves CSV correctly", {
  # Load dataset
  data("physeq_16SOTU", package = "DspikeIn")
  absolute_count <- microbiome::meta(physeq_16SOTU)
  
  # Check dataset validity
  expect_true(nrow(absolute_count) > 0, info = "Dataset should have rows.")
  expect_true(ncol(absolute_count) > 0, info = "Dataset should have columns.")
  
  # Create temporary CSV file path
  temp_csv <- tempfile(fileext = ".csv")
  
  # Run the function (writes CSV only)
  summary_table <- calculate_summary_stats_table(absolute_count, output_path = temp_csv)
  
  # --- Expected behavior ---
  expect_true(file.exists(temp_csv), info = "CSV file should be created.")
  
  # Verify CSV content
  summary_data <- read.csv(temp_csv)
  expect_true(nrow(summary_data) > 0, info = "CSV file should have rows.")
  expect_true(ncol(summary_data) > 0, info = "CSV file should have columns.")
  
  # Cleanup temporary file
  unlink(temp_csv)
})
