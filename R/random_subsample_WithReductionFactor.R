#' Random Subsampling with Reduction Factor
#'
#' This function performs random subsampling on the OTU table of a phyloseq object,
#' reducing the counts of each ASV (Amplicon Sequence Variant) by a specified reduction factor.
#' The resulting subsampled phyloseq object is saved to a specified output file.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param reduction_factor A numeric value specifying the factor by which to reduce the counts. Default is 3.
#' @param output_file A character string specifying the output file name for the subsampled phyloseq object. Default is "Less_subsampled_physeq.rds".
#' @return A phyloseq object with the subsampled OTU table.
#' @examples
#' # Perform random subsampling with a reduction factor of 10
#' red <- random_subsample_WithReductionFactor(spiked_16S, reduction_factor = 10)
#' # Summarize the subsampled phyloseq object
#' summ_phyloseq_sampleID(red)
#' @export
random_subsample_WithReductionFactor <- function(physeq, reduction_factor = 3, output_file = "Less_subsampled_physeq.rds") {
  require(phyloseq)
  
  # Get the OTU table
  otu_table_df <- as.data.frame(otu_table(physeq))
  
  # Loop through each sample
  for (i in 1:ncol(otu_table_df)) {
    # Get the counts for each ASV in the current sample
    asv_counts <- otu_table_df[, i]
    # Perform random subsampling for each ASV
    subsampled_counts <- numeric(length(asv_counts))
    for (j in seq_along(asv_counts)) {
      subsampled_counts[j] <- max(0, round(asv_counts[j] / reduction_factor))
    }
    
    # Update the OTU table with the subsampled counts
    otu_table_df[, i] <- subsampled_counts
  }
  
  # Convert the modified otu_table back to otu_table object
  otu_table_mod <- otu_table(otu_table_df, taxa_are_rows = TRUE)
  subsampled_physeq <- phyloseq(otu_table_mod, sample_data(physeq), tax_table(physeq))
  
  # Save the subsampled phyloseq object
  saveRDS(subsampled_physeq, file = output_file)
  cat("Less_Subsampled phyloseq object saved to:", output_file, "\n")
  
  return(subsampled_physeq)
}

# Example usage:
# Perform random subsampling with a reduction factor of 10
# red <- random_subsample_WithReductionFactor(spiked_16S, reduction_factor = 10)
# Summarize the subsampled phyloseq object
# summ_phyloseq_sampleID(red)
