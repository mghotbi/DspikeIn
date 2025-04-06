library(DspikeIn)
library(testthat)
library(phyloseq)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)

test_that("convert_to_absolute_counts works correctly with phyloseq and TSE", {
  output_file <- "absolute_counts.csv"

  if (file.exists(output_file)) file.remove(output_file)

  data("physeq_16SOTU", package = "DspikeIn")

  sample_names <- phyloseq::sample_names(physeq_16SOTU)
  scaling_factors <- setNames(rep(1000, length(sample_names)), sample_names)

  #  Run function on phyloseq object
  result_physeq <- convert_to_absolute_counts(physeq_16SOTU, scaling_factors)

  expect_type(result_physeq, "list")
  expect_true(all(c("absolute_counts", "obj_adj") %in% names(result_physeq)))

  expect_s3_class(result_physeq$absolute_counts, "data.frame")
  expect_s4_class(result_physeq$obj_adj, "phyloseq")

  #  Ensure absolute counts are non-negative
  expect_true(all(result_physeq$absolute_counts >= 0))

  expect_true(file.exists(output_file))

  #  Clean up after test
  file.remove(output_file)

  # =============================
  #  Test with TreeSummarizedExperiment (TSE)
  # =============================

  #  Convert phyloseq to TSE
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  sample_names_tse <- colnames(SummarizedExperiment::assay(tse_16SOTU))
  scaling_factors_tse <- setNames(rep(1000, length(sample_names_tse)), sample_names_tse)

  result_tse <- convert_to_absolute_counts(tse_16SOTU, scaling_factors_tse)

  #   result is a list
  expect_type(result_tse, "list")
  expect_true(all(c("absolute_counts", "obj_adj") %in% names(result_tse)))

  expect_s3_class(result_tse$absolute_counts, "data.frame")
  expect_s4_class(result_tse$obj_adj, "TreeSummarizedExperiment")

  #  Ensure absolute counts are non-negative
  expect_true(all(result_tse$absolute_counts >= 0))

  expect_true(file.exists(output_file))

  file.remove(output_file)

  # =============================
  #  Test Error Handling
  # =============================

  expect_error(
    convert_to_absolute_counts("not_a_physeq", scaling_factors),
    regexp = " Input must be a `phyloseq` or `TreeSummarizedExperiment` object."
  )

  #  Mismatched sample names in scaling factors
  invalid_scaling_factors <- setNames(rep(1000, length(sample_names)), paste0(sample_names, "_wrong"))
  expect_error(
    convert_to_absolute_counts(physeq_16SOTU, invalid_scaling_factors),
    regexp = " Error: Some sample names in scaling_factors do not match OTU table sample names."
  )
})
