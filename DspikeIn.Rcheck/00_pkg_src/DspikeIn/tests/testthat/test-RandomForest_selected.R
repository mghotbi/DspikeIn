# Load required packages for testing
library(testthat)
library(phyloseq)
library(randomForest)
library(DspikeIn)

test_that("RandomForest_selected works correctly with phyloseq", {
  data("physeq_16SOTU", package = "DspikeIn")

  result_phyloseq <- RandomForest_selected(
    physeq = physeq_16SOTU,
    response_var = "Host.genus",
    na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet")
  )

  expect_s4_class(result_phyloseq, "phyloseq")
  expect_gt(ntaxa(result_phyloseq), 0)

  # Clean up any auto-saved output
  unlink(list.files(pattern = "\\.csv$|\\.rds$", full.names = TRUE), force = TRUE)
})
