#' Proportionally Adjust Abundance
#'
#' This function normalizes the abundance data in a phyloseq object by adjusting each sample's counts
#' based on a total value, typically the maximum total sequence count across all samples. The adjusted
#' counts are then rounded to the nearest integer.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param output_file A character string specifying the output file name for the adjusted phyloseq object. Default is "proportion_adjusted_physeq.rds".
#' @return A phyloseq object with the proportionally adjusted and rounded abundance data.
#' @examples
#' # Proportionally adjust the abundance data
#' normalized_physeq <- proportion_adj(physeq, output_file = "proportion_adjusted_physeq.rds")
#' @export
proportion_adj <- function(physeq, output_file = "proportion_adjusted_physeq.rds") {
  # Normalize total sequence counts
  normf <- function(x, tot = max(sample_sums(physeq))) {
    tot * x / sum(x)
  }
  
  # Apply normalization to sample counts/abundance
  physeq <- transform_sample_counts(physeq, normf)
  
  # Round the abundance counts within the OTU table
  otu_table(physeq) <- round(otu_table(physeq), digits = 0)
  
  # Save the normalized and rounded phyloseq object
  saveRDS(physeq, file = output_file)
  cat("Normalized and rounded phyloseq object saved to:", output_file, "\n")
  
  # Return the normalized and rounded phyloseq obj
  return(physeq)
}

# Example usage:
# Proportionally adjust the abundance data
# normalized_physeq <- proportion_adj(physeq, output_file = "proportion_adjusted_physeq.rds")
