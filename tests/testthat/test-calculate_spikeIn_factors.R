testthat::skip_on_bioc()

library(testthat)
library(DspikeIn)

test_that("calculate_spike_percentage runs without errors", {
  data("physeq_16SOTU", package = "DspikeIn")

  merged_spiked_species <- c("Tetragenococcus_halophilus")

  result <- calculate_spike_percentage(
    obj = physeq_16SOTU,
    merged_spiked_species = merged_spiked_species,
    passed_range = c(0.1, 20)
  )

  expect_s3_class(result, "data.frame")
  expect_true("Percentage" %in% colnames(result))
  expect_true("Result" %in% colnames(result))

  # Clean up in case function writes to file
  unlink(list.files(pattern = "\\.csv$|\\.rda$", full.names = TRUE), force = TRUE)
})
