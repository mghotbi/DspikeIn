library(testthat)
library(DspikeIn)
library(phyloseq)
library(flextable)
library(dplyr)
library(utils)

test_that("calculate_spike_percentage() correctly computes percentages and saves files", {
  data("physeq_16SOTU", package = "DspikeIn")

  expect_true(ntaxa(physeq_16SOTU) > 0, "phyloseq object should have taxa.")
  expect_true(nsamples(physeq_16SOTU) > 0, "phyloseq object should have samples.")

  # Define spiked species for testing
  merged_spiked_species <- c("Tetragenococcus_halophilus")

  temp_docx <- tempfile(fileext = ".docx")
  temp_csv <- sub(".docx", ".csv", temp_docx)

  # Run the function with the corrected argument name
  result <- calculate_spike_percentage(
    physeq_16SOTU,
    merged_spiked_species = merged_spiked_species,
    output_file = temp_docx, # <- this was output_path before
    passed_range = c(0.1, 20)
  )

  expect_s3_class(result, "data.frame")

  expected_cols <- c("Sample", "Total_Reads", "Spiked_Reads", "Percentage", "Result")
  expect_true(all(expected_cols %in% colnames(result)), "Output should contain expected columns.")

  expect_type(result$Percentage, "double")

  expect_true(
    any(result$Result == "passed") | any(result$Result == "failed"),
    "Samples should be classified as 'passed' or 'failed'."
  )

  expect_true(file.exists(temp_docx), "Word document should be created.")
  expect_true(file.exists(temp_csv), "CSV file should be created.")

  # Read the CSV file and check structure
  saved_data <- read.csv(temp_csv)
  expect_true(nrow(saved_data) > 0, "CSV file should have rows.")
  expect_true(ncol(saved_data) > 0, "CSV file should have columns.")

  # Cleanup temporary files
  unlink(temp_docx)
  unlink(temp_csv)
})
