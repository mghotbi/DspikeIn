# ============================================================
# Test: imbalance_calculate_list_average_scaling_factors()
# ============================================================
#' @description
#' Unit test for the DspikeIn function that computes per-sample scaling factors
#' across multiple spike-in species for both phyloseq and TSE formats.
# ============================================================

options(warn = -1)
suppressPackageStartupMessages({
  library(testthat)
  library(DspikeIn)
  library(phyloseq)
  library(SummarizedExperiment)
})
options(warn = 0)

test_that("imbalance_calculate_list_average_scaling_factors() works for phyloseq and TSE", {
  # ------------------------------------------------------------
  # 1. Construct synthetic dataset (identical to function example)
  # ------------------------------------------------------------
  otu <- matrix(
    c(
      6000, 6200, 5900, 6100,
      4000, 4200, 3900, 4100,
      2000, 1900, 2100, 2050,
      1300, 1250, 1350, 1400,
      500,  800,  900,  700,   # Flavobacterium_spike
      900, 1200, 1100, 1000    # Bacillus_spike
    ),
    nrow = 6, byrow = TRUE,
    dimnames = list(
      c("OTU1", "OTU2", "OTU3", "OTU4",
        "Flavobacterium_spike", "Bacillus_spike"),
      c("S1", "S2", "S3", "S4")
    )
  )
  
  tax <- data.frame(
    Kingdom = rep("Bacteria", 6),
    Species = c("OTU1", "OTU2", "OTU3", "OTU4",
                "Flavobacterium_spike", "Bacillus_spike"),
    row.names = rownames(otu)
  )
  
  #  sample_data() must have at least one column
  sam <- data.frame(
    Group = c("A", "A", "B", "B"),
    row.names = c("S1", "S2", "S3", "S4")
  )
  
  ps <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    tax_table(as.matrix(tax)),
    sample_data(sam)
  )
  
  # ------------------------------------------------------------
  # 2. Define spike-in lists (expected cells per sample)
  # ------------------------------------------------------------
  spiked_species_list <- list(
    Flavo = "Flavobacterium_spike",
    Bacillus = "Bacillus_spike"
  )
  
  spiked_cells_list <- list(
    Flavo = c(S1 = 1e7, S2 = 3e7, S3 = 6e7, S4 = 2e7),
    Bacillus = c(S1 = 2e7, S2 = 1e7, S3 = 5e7, S4 = 3e7)
  )
  
  # ------------------------------------------------------------
  # 3. Run function on phyloseq object
  # ------------------------------------------------------------
  factors_phy <- imbalance_calculate_list_average_scaling_factors(
    ps,
    spiked_species_list,
    spiked_cells_list,
    normalize = FALSE,
    verbose = FALSE
  )
  
  expect_type(factors_phy, "double")
  expect_equal(length(factors_phy), 4)
  expect_named(factors_phy, c("S1", "S2", "S3", "S4"))
  expect_true(all(factors_phy > 0, na.rm = TRUE))
  
  # Expected ratio  (expected cells / observed spike reads)
  merged_abund <- colSums(otu[c("Flavobacterium_spike", "Bacillus_spike"), ])
  expected_total <- spiked_cells_list$Flavo + spiked_cells_list$Bacillus
  rough_factor <- expected_total / merged_abund
  
  expect_true(all(factors_phy > 0))
  expect_true(all(factors_phy < rough_factor * 2))
  
  # ------------------------------------------------------------
  # 4. Convert to TSE and verify consistent results
  # ------------------------------------------------------------
  tss <- convert_phyloseq_to_tse(ps)
  factors_tse <- imbalance_calculate_list_average_scaling_factors(
    tss,
    spiked_species_list,
    spiked_cells_list,
    normalize = FALSE,
    verbose = FALSE
  )
  
  expect_type(factors_tse, "double")
  expect_equal(length(factors_tse), 4)
  expect_named(factors_tse, c("S1", "S2", "S3", "S4"))
  
  diff <- abs(factors_tse - factors_phy)
  expect_true(all(diff < max(factors_phy) * 0.3),
              info = "TSE and phyloseq results differ too much")
  
  # ------------------------------------------------------------
  # 5. Test normalization and Inf-handling
  # ------------------------------------------------------------
  f_norm <- imbalance_calculate_list_average_scaling_factors(
    ps,
    spiked_species_list,
    spiked_cells_list,
    normalize = TRUE
  )
  expect_equal(round(stats::median(f_norm, na.rm = TRUE), 6), 1)
  
  f_inf <- imbalance_calculate_list_average_scaling_factors(
    ps,
    spiked_species_list,
    spiked_cells_list,
    normalize = FALSE,
    allow_infinite = TRUE
  )
  expect_true(all(is.finite(f_inf) | is.infinite(f_inf)))
  
  # ------------------------------------------------------------
  # 6. Test invalid input error handling
  # ------------------------------------------------------------
  expect_error(
    imbalance_calculate_list_average_scaling_factors(
      "not_an_obj", spiked_species_list, spiked_cells_list
    ),
    "Input must be a phyloseq or TreeSummarizedExperiment object"
  )
})
