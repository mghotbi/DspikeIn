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

<<<<<<< HEAD
  # ─── Test with phyloseq ──────────────────────────────────────
=======
  #  Test with phyloseq 
>>>>>>> bioc-rev1
  result_physeq <- convert_to_absolute_counts(physeq_16SOTU, scaling_factors)

  expect_type(result_physeq, "list")
  expect_named(result_physeq, c("absolute_counts", "obj_adj"))

  expect_s3_class(result_physeq$absolute_counts, "data.frame")
  expect_s4_class(result_physeq$obj_adj, "phyloseq")
  expect_true(all(result_physeq$absolute_counts >= 0))
  expect_true(file.exists(output_file))
  file.remove(output_file)

<<<<<<< HEAD
  # ─── Test with TreeSummarizedExperiment ──────────────────────
=======
  #  Test with TreeSummarizedExperiment 
>>>>>>> bioc-rev1
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  sample_names_tse <- colnames(SummarizedExperiment::assay(tse_16SOTU))
  scaling_factors_tse <- setNames(rep(1000, length(sample_names_tse)), sample_names_tse)

  result_tse <- convert_to_absolute_counts(tse_16SOTU, scaling_factors_tse)

  expect_type(result_tse, "list")
  expect_named(result_tse, c("absolute_counts", "obj_adj"))

  expect_s3_class(result_tse$absolute_counts, "data.frame")
  expect_s4_class(result_tse$obj_adj, "TreeSummarizedExperiment")
  expect_true(all(result_tse$absolute_counts >= 0))
  expect_true(file.exists(output_file))
  file.remove(output_file)

<<<<<<< HEAD
  # ─── Error handling: invalid object ──────────────────────────
=======
  #  Error handling: invalid object 
>>>>>>> bioc-rev1
  expect_error(
    convert_to_absolute_counts("not_a_physeq", scaling_factors),
    regexp = "must be a `phyloseq` or `TreeSummarizedExperiment` object"
  )

<<<<<<< HEAD
  # ─── Error handling: mismatched scaling factors ──────────────
=======
  #  Error handling: mismatched scaling factors 
>>>>>>> bioc-rev1
  bad_factors <- setNames(rep(1000, length(sample_names)), paste0(sample_names, "_x"))
  expect_error(
    convert_to_absolute_counts(physeq_16SOTU, bad_factors),
    regexp = "sample names in scaling_factors do not match OTU table sample names"
  )
})
