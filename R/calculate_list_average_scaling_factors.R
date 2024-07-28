#' Calculate Average Scaling Factors for Multiple Spiked Species with Equal or Different Spiked Cell Counts
#'
#' This function calculates scaling factors for multiple spiked species in a phyloseq object,
#' merges ASVs for each species if necessary, averages the scaling factors, and returns the
#' averaged scaling factors. Different spiked cell counts can be provided for each set of spiked species.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param spiked_species_list A list of character vectors, where each vector contains the spiked species to merge for each scaling factor calculation.
#' @param spiked_cells_list A numeric vector specifying the number of spiked cells for each set of spiked species.
#' @param merge_method The method to use for merging ASVs: "sum" or "max".
#' @return A numeric vector of averaged scaling factors.
#' @examples
#' # Example usage:
#' # Assume `physeq` is a phyloseq object and `spiked_species_list` is a list of character vectors.
#' # spiked_species_list <- list(c("Species1"), c("Species2"), c("Species3"), c("Species4"))
#' # spiked_cells_list <- c(1000, 1500, 2000, 2500)
#' # average_scaling_factors <- calculate_average_scaling_factors(physeq, spiked_species_list, spiked_cells_list, "sum")
#' @export
calculate_list_average_scaling_factors <- function(physeq, spiked_species_list, spiked_cells_list, merge_method = c("sum", "max")) {
  merge_method <- match.arg(merge_method)
  
  # Check if lengths of spiked_species_list and spiked_cells_list match
  if (length(spiked_species_list) != length(spiked_cells_list)) {
    stop("The length of spiked_species_list must match the length of spiked_cells_list")
  }
  
  # Initialize a list to store individual scaling factors
  scaling_factors_list <- list()
  
  # Loop through each set of spiked species and corresponding spiked cell count
  for (i in seq_along(spiked_species_list)) {
    # Pre-process species to merge ASVs
    processed_physeq <- Pre_processing_species(physeq, spiked_species_list[[i]], merge_method)
    
    # Calculate scaling factors for the processed phyloseq object
    result <- calculate_spikeIn_factors(processed_physeq, spiked_cells_list[i], spiked_species_list[[i]])
    scaling_factors_list[[length(scaling_factors_list) + 1]] <- result$scaling_factors
  }
  
  # Calculate the average of the scaling factors
  average_scaling_factors <- rowMeans(do.call(cbind, scaling_factors_list), na.rm = TRUE)
  
  return(average_scaling_factors)
}

# Example usage:
# Define the spiked species list for four different spiked species
# spiked_species_list <- list(
#   c("Methylobacterium_phyllostachyos"),
#   c("Methylorubrum_salsuginis"),
#   c("Bosea_massiliensis"),
#   c("Bacillus_decolorationis")
# )
# 
# # Example usage;
# # # Define the corresponding number of spiked cells for each spiked species 
# spiked_cells_list <- c(20000, 20000, 20000, 20000)
# 
# # Calculate the average scaling factors with ASV merging using the "sum" method
# average_scaling_factors <- calculate_list_average_scaling_factors(ps, spiked_species_list, spiked_cells_list, "sum")
# 
# # Convert the relative counts to absolute counts
# absolute <- convert_to_absolute_counts(ps, average_scaling_factors)
# absolute_counts <- absolute$absolute_counts
# physeq_absolute_abundance <- absolute$physeq_obj

