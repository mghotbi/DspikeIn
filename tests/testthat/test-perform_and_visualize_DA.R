library(testthat)
library(phyloseq)
library(DspikeIn)
library(ggplot2)
library(TreeSummarizedExperiment)

test_that("perform_and_visualize_DA returns expected outputs", {
  data("physeq_16SOTU", package = "DspikeIn")

  #  Define test parameters
  method <- "DESeq2"
  group_var <- "Diet"
  contrast <- c("Insectivore", "Carnivore")
  significance_level <- 0.05

  result_phyloseq <- tryCatch(
    {
      suppressWarnings(perform_and_visualize_DA(
        obj = physeq_16SOTU,
        method = method,
        group_var = group_var,
        contrast = contrast,
        significance_level = significance_level
      ))
    },
    error = function(e) {
      cat("\n ERROR in perform_and_visualize_DA: ", conditionMessage(e), "\n")
      NULL
    }
  )

  #  1. Ensure function did not return NULL
  expect_false(is.null(result_phyloseq))

  #  2. Ensure function returns a list
  expect_type(result_phyloseq, "list")

  #  3. Ensure expected keys exist in output
  expect_true(all(c("results", "obj_significant", "plot") %in% names(result_phyloseq)))

  #  4. Ensure `results` is a dataframe and has correct structure
  expect_s3_class(result_phyloseq$results, "data.frame")
  expect_true(all(c("logFC", "pvalue", "FDR") %in% colnames(result_phyloseq$results)))

  #  5. Ensure `plot` is a ggplot object
  expect_s3_class(result_phyloseq$plot, "ggplot")

  #  6. Ensure `obj_significant` is valid (phyloseq or TSE)
  expect_true(
    inherits(result_phyloseq$obj_significant, "phyloseq") ||
      inherits(result_phyloseq$obj_significant, "TreeSummarizedExperiment")
  )

  #  7. Repeat test for TreeSummarizedExperiment input
  tse_test <- DspikeIn::convert_phyloseq_to_tse(physeq_16SOTU)

  result_TSE <- suppressWarnings(perform_and_visualize_DA(
    obj = tse_test,
    method = method,
    group_var = group_var,
    contrast = contrast,
    significance_level = significance_level
  ))

  #  Ensure output remains valid for TSE input
  expect_type(result_TSE, "list")
  expect_true(inherits(result_TSE$obj_significant, "TreeSummarizedExperiment"))
})
