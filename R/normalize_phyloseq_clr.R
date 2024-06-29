#' Centered Log-Ratio (CLR)
#'
#' This function performs CLR normalization on a phyloseq object using log-ratio transformation.
#' It handles zero, negative, and NA values by removing them and informing the user of the modifications.
#'
#' @param ps A phyloseq object.
#' @param feature_category Character string specifying which features to use as the divisor for the geometric mean calculation.
#' @param min_counts Minimum number of counts required for a sample to be retained.
#' @return A phyloseq object with CLR normalization.
#' @export
#' @examples
#' data(GlobalPatterns)
#' ps_clr <- normalize_phyloseq_clr(GlobalPatterns, feature_category = "iqlr", min_counts = 1000)
normalize_phyloseq_clr <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha"), min_counts = 1) {
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
  
  gm_mean <- function(x) exp(sum(log(x[x > 0])) / length(x))
  phyloseq::transform_sample_counts(ps, function(x) log(x / gm_mean(x)))
}

# Example;
# physeq <- phy
# physeq_clr <- normalize_phyloseq_clr(physeq, feature_category = "zero", min_counts = 1000)
