testthat::skip_on_bioc()

library(testthat)
library(phyloseq)
library(TreeSummarizedExperiment)
library(DspikeIn)

test_that("remove_zero_negative_count_samples removes zero/negative samples and adds pseudocounts", {
  # Load test data
  data("physeq_16SOTU", package = "DspikeIn")
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Introduce zero/negative counts in a sample
  otu_physeq <- get_otu_table(physeq_16SOTU)
  otu_physeq[, 1] <- 0 # Set the first sample to zero
  physeq_16SOTU <- phyloseq::phyloseq(
    otu_table(otu_physeq, taxa_are_rows = TRUE),
    sample_data(physeq_16SOTU),
    tax_table(physeq_16SOTU)
  )

  # Process with function
  cleaned_physeq <- remove_zero_negative_count_samples(physeq_16SOTU, pseudocount = 1)

  expect_true(inherits(cleaned_physeq, "phyloseq"))

  # Ensure zero-count samples are removed
  expect_false(any(phyloseq::sample_sums(cleaned_physeq) == 0))

  # Ensure pseudocount was added
  expect_true(min(get_otu_table(cleaned_physeq)) > 0)

  # Repeat for TSE
  otu_tse <- get_otu_table(tse_16SOTU)
  otu_tse[, 1] <- 0 # Set the first sample to zero
  SummarizedExperiment::assay(tse_16SOTU) <- otu_tse

  cleaned_tse <- remove_zero_negative_count_samples(tse_16SOTU, pseudocount = 1)

  # Ensure output is still a TSE object
  expect_true(inherits(cleaned_tse, "TreeSummarizedExperiment"))

  # Ensure zero-count samples are removed
  expect_false(any(colSums(get_otu_table(cleaned_tse)) == 0))

  # Ensure pseudocount was added
  expect_true(min(get_otu_table(cleaned_tse)) > 0)
})
