library(testthat)
library(phyloseq)
library(TreeSummarizedExperiment)
library(DspikeIn)

test_that("Pre_processing_species correctly processes a phyloseq object", {
  # Load test data
  data("physeq_16SOTU", package = "DspikeIn")

  # species to merge
  species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")

  # Run function
  merged_physeq <- Pre_processing_species(physeq_16SOTU, species_name, merge_method = "sum")

  # Ensure output is a phyloseq object
  expect_true(inherits(merged_physeq, "phyloseq"))

  expect_true(!is.null(otu_table(merged_physeq)))
  expect_true(!is.null(tax_table(merged_physeq)))
  expect_true(!is.null(sample_data(merged_physeq)))

  # Ensure phylogenetic tree and reference sequences remain
  expect_true(!is.null(phy_tree(merged_physeq)))
  expect_true(!is.null(refseq(merged_physeq)))

  # 🔻 Clean up any files generated (csv, rda)
  unlink(list.files(pattern = "\\.csv$|\\.rda$", full.names = TRUE), force = TRUE)
})
