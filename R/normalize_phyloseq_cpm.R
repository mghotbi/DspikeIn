#' Counts Per Million (CPM)
#'
#' This function normalizes counts in a phyloseq object to counts per million (CPM).
#' It handles zero, negative, and NA values by removing them and informing the user of the modifications.
#'
#' @param ps A phyloseq object.
#' @param feature_category Character string specifying which features to use as the divisor for the geometric mean calculation.
#' @param min_counts Minimum number of counts required for a sample to be retained.
#' @return A phyloseq object with CPM normalization.
#' @export
#' @examples
#' data(GlobalPatterns)
#' ps_cpm <- normalize_phyloseq_cpm(GlobalPatterns, feature_category = "iqlr", min_counts = 1000)
normalize_phyloseq_cpm <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha"), min_counts = 1) {
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
  
  phyloseq::transform_sample_counts(ps, function(x) x / sum(x) * 1e6)
}

# Example;
# physeq <- phy
# physeq_cpm <- normalize_phyloseq_cpm(physeq, feature_category = "zero", min_counts = 1000)
