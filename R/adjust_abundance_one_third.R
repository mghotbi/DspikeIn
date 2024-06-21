#' Adjust Abundance by One-Third
#'
#' This function normalizes the abundance data in a phyloseq object by dividing each value by a specified factor.
#' The normalization process is based on the volume of DNA processed for sequencing, such as 16S rRNA amplicons.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param factor A numeric value specifying the factor by which to divide the abundance data. Default is 3.
#' @param output_file A character string specifying the output file name for the adjusted phyloseq object. Default is "physeq_adj_scaled.rds".
#' @return A phyloseq object with the adjusted abundance data.
#' @examples
#' # Adjust the abundance data by dividing each value by 3
#' adjusted_physeq <- adjust_abundance_one_third(physeq16S, factor = 3)
#' @export
adjust_abundance_one_third <- function(physeq, factor = 3, output_file = "physeq_adj_scaled.rds") {
  print("Starting adjustment process...")
  
  # Check if the OTU table is already a matrix
  if (class(physeq@otu_table) == "matrix") {
    print("OTU table is already a matrix. Performing division...")
    # Divide the OTU table by the specified factor
    physeq@otu_table <- physeq@otu_table / factor
  } else {
    print("Converting OTU table to matrix...")
    # Convert the OTU table to a matrix and divide by the specified factor
    physeq@otu_table <- as.matrix(physeq@otu_table) / factor
  }
  
  print("Saving modified phyloseq object...")
  # Save the adjusted phyloseq object to a file
  saveRDS(physeq, file = output_file)
  cat("Modified phyloseq object saved as:", output_file, "\n")
  
  print("Adjustment complete.")
  
  return(physeq)
}

# Example usage:
# Adjust the abundance data by dividing each value by 3
# adjusted_physeq <- adjust_abundance_one_third(physeq16S, factor = 3)
