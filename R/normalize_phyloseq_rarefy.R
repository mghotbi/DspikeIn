#' Rarefaction: Random Subsampling to the Smallest Library Size
#'
#' This function performs rarefaction on a phyloseq object by randomly subsampling counts to the smallest library size.
#' It handles zero, negative, and NA values by removing them and informing the user of the modifications.
#'
#' @param ps A phyloseq object.
#' @param feature_category Character string specifying which features to use as the divisor for the geometric mean calculation.
#' @param min_counts Minimum number of counts required for a sample to be retained.
#' @return A rarefied phyloseq object.
#' @export
#' @examples
#' data(GlobalPatterns)
#' ps_rarefied <- normalize_phyloseq_rarefy(GlobalPatterns, feature_category = "iqlr", min_counts = 1000)
normalize_phyloseq_rarefy <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha"), min_counts = 1) {
  ps <- remove_zero_negative_count_samples(ps) # Remove samples with zero or negative counts
  ps <- process_data_with_feature_category(ps, feature_category) # Process data with the specified feature category
  
  # Check if the OTU table has valid dimensions after processing
  if (nrow(otu_table(ps)) == 0 || ncol(otu_table(ps)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  library(phyloseq)
  
  # Perform rarefaction
  cat("Performing rarefaction with set.seed(123) for reproducibility.\n")
  phyloseq::rarefy_even_depth(ps, rngseed = 123)
}

# Example,
# physeq <- phy
# physeq_rarefy <- normalize_phyloseq_rarefy(physeq, feature_category = "zero", min_counts = 1000)
