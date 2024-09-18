#' Adjust Abundance by One-Third
#'
#' This function normalizes the abundance data in a phyloseq object by dividing each value by a specified factor.
#' The normalization process is based on the volume of DNA processed for sequencing, such as 16S rRNA amplicons.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param factor A numeric value specifying the factor by which to divide the abundance data. Default is 3.
#' @param output_file A character string specifying the output file name for the adjusted phyloseq object. Default is NULL.
#' @return A phyloseq object with the adjusted abundance data.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Adjust the abundance data by dividing each value by 3
#'   adjusted_physeq <- adjust_abundance_one_third(physeq16S, factor = 3)
#' }
#' }
#' @export
adjust_abundance_one_third <- function(physeq, factor = 3, output_file = NULL) {
  suppressMessages({
    message("Starting adjustment process...")
    
    # Check if the OTU table is already a matrix
    otu_table_matrix <- phyloseq::otu_table(physeq)
    
    if (is.matrix(otu_table_matrix)) {
      message("OTU table is already a matrix. Performing division...")
      # Divide the OTU table by the specified factor
      otu_table_matrix <- otu_table_matrix / factor
    } else {
      message("Converting OTU table to matrix...")
      # Convert the OTU table to a matrix and divide by the specified factor
      otu_table_matrix <- as.matrix(otu_table_matrix) / factor
    }
    
    # Update the OTU table in the phyloseq object
    phyloseq::otu_table(physeq) <- otu_table_matrix
    
    # Save the adjusted phyloseq object to an output file if specified
    if (!is.null(output_file)) {
      message("Saving modified phyloseq object to: ", output_file)
      saveRDS(physeq, file = output_file)
    }
    
    return(physeq)
  })
}
# Example usage:
# Adjust the abundance data by dividing each value by 3
# adjusted_physeq <- adjust_abundance_one_third(physeq_16SASV, factor = 3)
