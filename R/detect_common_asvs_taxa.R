#' Detect Common ASVs and Taxa from Multiple Phyloseq Objects
#'
#' This function identifies the common ASVs (Amplicon Sequence Variants) and common taxa detected across multiple phyloseq objects.
#' It extracts ASVs and taxa from the provided phyloseq objects, finds the common ones, and optionally saves the results to CSV and RDS files.
#'
#' @param phyloseq_list A list of phyloseq objects.
#' @param output_common_asvs_csv A character string specifying the path to save the common ASVs as a CSV file. Default is "common_asvs.csv".
#' @param output_common_asvs_rds A character string specifying the path to save the common ASVs as an RDS file. Default is "common_asvs.rds".
#' @param output_common_taxa_csv A character string specifying the path to save the common taxa as a CSV file. Default is "common_taxa.csv".
#' @param output_common_taxa_rds A character string specifying the path to save the common taxa as an RDS file. Default is "common_taxa.rds".
#' @param return_as_df A logical indicating whether to return the results as data frames. Default is FALSE, meaning the results are returned as phyloseq objects.
#' @return A list containing the phyloseq objects or data frames of common ASVs and common taxa.
#' @examples
#' # Example usage:
#' # results <- detect_common_asvs_taxa(list(physeq1, physeq2, physeq3), return_as_df = TRUE)
#' # common_asvs_df <- results$common_asvs
#' # common_taxa_df <- results$common_taxa
#' @export
detect_common_asvs_taxa <- function(phyloseq_list, 
                                    output_common_asvs_csv = "common_asvs.csv", 
                                    output_common_asvs_rds = "common_asvs.rds", 
                                    output_common_taxa_csv = "common_taxa.csv", 
                                    output_common_taxa_rds = "common_taxa.rds",
                                    return_as_df = FALSE) {
  # Check if the list has at least two phyloseq objects
  if (length(phyloseq_list) < 2) {
    stop("At least two phyloseq objects are required.")
  }
  
  # Check if tax_table slot is not empty for each phyloseq object
  for (i in seq_along(phyloseq_list)) {
    if (is.null(tax_table(phyloseq_list[[i]]))) {
      stop(paste("tax_table slot is empty in phyloseq object at position", i))
    }
  }
  
  # Extract ASVs and taxa from all phyloseq objects
  asv_lists <- lapply(phyloseq_list, function(x) {
    rownames(otu_table(x))
  })
  
  taxa_lists <- lapply(phyloseq_list, function(x) {
    rownames(tax_table(x))
  })
  
  # Find common ASVs and taxa
  common_asvs <- Reduce(intersect, asv_lists)
  common_taxa <- Reduce(intersect, taxa_lists)
  
  # Initialize result list
  result <- list(common_asvs = NULL, common_taxa = NULL)
  
  # Create data frames or phyloseq objects for common ASVs if any are found
  if (length(common_asvs) > 0) {
    common_asvs_phyloseq <- prune_taxa(common_asvs, phyloseq_list[[1]])
    common_asvs_df <- psmelt(common_asvs_phyloseq)
    
    write.csv(common_asvs_df, output_common_asvs_csv, row.names = FALSE)
    cat("Common ASVs saved to:", output_common_asvs_csv, "\n")
    
    saveRDS(common_asvs_phyloseq, file = output_common_asvs_rds)
    cat("Common ASVs saved to:", output_common_asvs_rds, "\n")
    
    result$common_asvs <- if (return_as_df) common_asvs_df else common_asvs_phyloseq
  } else {
    cat("No common ASVs found.\n")
  }
  
  # Create data frames or phyloseq objects for common taxa if any are found
  if (length(common_taxa) > 0) {
    common_taxa_phyloseq <- prune_taxa(common_taxa, phyloseq_list[[1]])
    common_taxa_df <- psmelt(common_taxa_phyloseq)
    
    write.csv(common_taxa_df, output_common_taxa_csv, row.names = FALSE)
    cat("Common taxa saved to:", output_common_taxa_csv, "\n")
    
    saveRDS(common_taxa_phyloseq, file = output_common_taxa_rds)
    cat("Common taxa saved to:", output_common_taxa_rds, "\n")
    
    result$common_taxa <- if (return_as_df) common_taxa_df else common_taxa_phyloseq
  } else {
    cat("No common taxa found.\n")
  }
  
  return(result)
}


# Example usage:
#results <- detect_common_asvs_taxa(  list(rf_physeq, FT, core),  return_as_df = TRUE)

# Access the results
#common_asvs <- results$common_asvs
#common_taxa <- results$common_taxa

