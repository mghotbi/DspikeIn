#' Relative Log Expression (RLE)
#'
#' This function performs RLE normalization on a phyloseq object using a pseudo-reference sample.
#' It handles zero, negative, and NA values by removing them and informing the user of the modifications.
#'
#' @param ps A phyloseq object.
#' @param feature_category Character string specifying which features to use as the divisor for the geometric mean calculation.
#' @param min_counts Minimum number of counts required for a sample to be retained.
#' @return A phyloseq object with RLE normalization.
#' @export
#' @examples
#' data(GlobalPatterns)
#' ps_rle <- normalize_phyloseq_rle(GlobalPatterns, feature_category = "iqlr", min_counts = 1000)
normalize_phyloseq_rle <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha"), min_counts = 1) {
  cat("Initial OTU table dimensions: ", dim(otu_table(ps)), "\n")
  
  ps <- remove_zero_negative_count_samples(ps)
  cat("OTU table dimensions after removing zero/negative count samples: ", dim(otu_table(ps)), "\n")
  
  ps <- process_data_with_feature_category(ps, feature_category)
  otu_table_dim <- dim(otu_table(ps))
  cat("OTU table dimensions after processing feature category '", feature_category, "': ", otu_table_dim, "\n")
  
  if (otu_table_dim[1] == 0 || otu_table_dim[2] == 0) {
    stop("OTU table has zero dimensions after processing. Ensure the data contains valid counts.")
  }
  
  counts <- phyloseq::otu_table(ps)
  
  # Filter out taxa with all zero counts
  non_zero_taxa <- rowSums(counts) > 0
  counts <- counts[non_zero_taxa, ]
  
  cat("Summary of counts after filtering out zero-count taxa:\n")
  print(summary(counts))
  
  if (any(rowSums(counts) == 0)) {
    stop("OTU table contains samples with zero counts after processing. Ensure the data contains valid counts.")
  }
  
  # Add a small pseudocount to avoid issues with zero counts
  counts <- counts + 1
  
  # Normalize using RLE
  norm_counts <- norm_rle(counts)
  
  cat("Normalized counts dimensions: ", dim(norm_counts), "\n")
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(norm_counts, taxa_are_rows = TRUE)
  
  cat("Final OTU table dimensions after RLE normalization: ", dim(otu_table(ps)), "\n")
  
  return(ps)
}

# Example;
# physeq <- phy
# physeq_rle <- normalize_phyloseq_rle(physeq, feature_category = "zero", min_counts = 1000)
