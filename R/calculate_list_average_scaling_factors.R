#' Calculate Average Scaling Factors for Multiple Spiked Species
#'
#' This function calculates scaling factors for multiple spiked species in a \code{phyloseq} object.
#' It merges ASVs/OTUs for each species if necessary, averages the scaling factors, and returns
#' the averaged scaling factors for each OTU. Different spiked cell counts can be provided for each set of spiked species.
#' If an OTU is not associated with any spiked species, a default scaling factor (1) is assigned.
#' Scaling factors are rounded to the specified number of decimal places.
#'
#' @param physeq A \code{phyloseq} object containing microbial abundance data.
#' @param spiked_species_list A list of character vectors. Each vector contains the spiked species
#' (by taxon names) to be merged for calculating scaling factors for each group.
#' @param spiked_cells_list A numeric vector specifying the number of spiked cells corresponding to each group in `spiked_species_list`.
#' @param merge_method A character string specifying how to merge ASVs/OTUs for each group of spiked species.
#' Accepted values are \code{"sum"} or \code{"max"}. Default is \code{"sum"}.
#' @return A numeric vector of averaged and rounded scaling factors for each OTU in the \code{phyloseq} object.
#' OTUs not associated with any spiked species are assigned a default scaling factor of 1.
#' @examples
#' \dontrun{
#' # Example usage:
#' spiked_species_list <- list(
#'   c("Pseudomonas aeruginosa"),
#'   c("Escherichia coli"),
#'   c("Clostridium difficile")
#' )
#' spiked_cells_list <- c(10000, 20000, 15000)
#' scaling_factors <- calculate_list_average_scaling_factors(physeq, spiked_species_list, spiked_cells_list, merge_method = "sum")
#' print(scaling_factors)
#' }
#' @importFrom phyloseq taxa_names otu_table
#' @export
calculate_list_average_scaling_factors <- function(physeq, spiked_species_list, spiked_cells_list, merge_method = c("sum", "max")) {
  
  # Ensure correct merge_method input
  merge_method <- match.arg(merge_method)
  
  # Suppress messages for loading required packages
  suppressMessages({
    requireNamespace("phyloseq", quietly = TRUE)
  })
  
  # Check if lengths of spiked_species_list and spiked_cells_list match
  if (length(spiked_species_list) != length(spiked_cells_list)) {
    stop("The length of spiked_species_list must match the length of spiked_cells_list.")
  }
  
  # Initialize vectors for combined scaling factors and counts
  otu_names <- phyloseq::taxa_names(physeq)  # Get OTU names
  scaling_factors_total <- rep(0, length(otu_names))  # Initialize with zeros
  count_contributions <- rep(0, length(otu_names))  # To track scaling factor contributions
  
  # Loop through each group of spiked species and calculate scaling factors
  for (i in seq_along(spiked_species_list)) {
    
    # Pre-process species to merge ASVs/OTUs if necessary
    processed_physeq <- Pre_processing_species_list(physeq, spiked_species_list[[i]], merge_method)
    
    # Calculate the total observed abundance for the spiked species in the phyloseq object
    otu_table_data <- as(phyloseq::otu_table(processed_physeq), "matrix")
    total_abundance_spiked <- sum(otu_table_data)
    
    # Calculate the scaling factor based on the ratio of spiked cells to the observed total abundance
    scaling_factor <- spiked_cells_list[i] / total_abundance_spiked
    
    # Get the OTU names for the spiked species
    spiked_otu_names <- phyloseq::taxa_names(processed_physeq)
    
    # Find indices of the OTUs in spiked species set
    matching_otu_indices <- match(spiked_otu_names, otu_names)
    
    # Remove NA values from matching_otu_indices
    valid_indices <- !is.na(matching_otu_indices)
    if (any(valid_indices)) {
      # Accumulate the scaling factors into the total scaling factors vector
      scaling_factors_total[matching_otu_indices[valid_indices]] <- scaling_factors_total[matching_otu_indices[valid_indices]] + scaling_factor
      count_contributions[matching_otu_indices[valid_indices]] <- count_contributions[matching_otu_indices[valid_indices]] + 1
    }
  }
  
  # Average the scaling factors (divide by number of contributions per OTU)
  average_scaling_factors <- scaling_factors_total / count_contributions
  
  # Replace NAs with default scaling factor of 1 for OTUs with no contributions
  average_scaling_factors[count_contributions == 0] <- NA
  average_scaling_factors[is.na(average_scaling_factors)] <- 1
  
  # Round the scaling factors to 2 decimal places
  average_scaling_factors <- round(average_scaling_factors, digits = 2)
  
  return(average_scaling_factors)
}


# Example usage:
# # Step 1: Create a taxonomy table 
# taxa_data <- data.frame(
#   OTUID = c("ASV1", "ASV2", "ASV3",    
#             "ASV4", "ASV5", "ASV6",   
#             "ASV7", "ASV8", "ASV9",    
#             "ASV10", "ASV11", "ASV12", 
#             "ASV13", "ASV14", "ASV15", 
#             "ASV16", "ASV17", "ASV18"),
#   Kingdom = rep("Bacteria", 18),
#   Phylum = c("Proteobacteria", "Proteobacteria", "Proteobacteria", 
#              "Proteobacteria", "Proteobacteria", "Proteobacteria", 
#              "Firmicutes", "Firmicutes", "Firmicutes", 
#              "Proteobacteria", "Firmicutes", "Firmicutes", 
#              "Firmicutes", "Firmicutes", "Firmicutes", 
#              "Bacteroidota", "Bacteroidota", "Proteobacteria"),
#   Class = c("Gammaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria", 
#             "Gammaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria", 
#             "Clostridia", "Clostridia", "Clostridia", 
#             "Gammaproteobacteria", "Bacilli", "Bacilli", 
#             "Bacilli", "Bacilli", "Clostridia", 
#             "Bacteroidia", "Bacteroidia", "Epsilonproteobacteria"),
#   Order = c("Pseudomonadales", "Pseudomonadales", "Pseudomonadales", 
#             "Enterobacterales", "Enterobacterales", "Enterobacterales", 
#             "Clostridiales", "Clostridiales", "Clostridiales", 
#             "Enterobacterales", "Lactobacillales", "Bacillales", 
#             "Bacillales", "Bacillales", "Clostridiales", 
#             "Bacteroidales", "Bacteroidales", "Campylobacterales"),
#   Family = c("Pseudomonadaceae", "Pseudomonadaceae", "Pseudomonadaceae", 
#              "Enterobacteriaceae", "Enterobacteriaceae", "Enterobacteriaceae", 
#              "Clostridiaceae", "Clostridiaceae", "Clostridiaceae", 
#              "Enterobacteriaceae", "Enterococcaceae", "Staphylococcaceae", 
#              "Listeriaceae", "Bacillaceae", "Clostridiaceae", 
#              "Bacteroidaceae", "Bacteroidaceae", "Helicobacteraceae"),
#   Genus = c("Pseudomonas", "Pseudomonas", "Pseudomonas", 
#             "Escherichia", "Escherichia", "Escherichia", 
#             "Clostridium", "Clostridium", "Clostridium", 
#             "Salmonella", "Enterococcus", "Staphylococcus", 
#             "Listeria", "Bacillus", "Lactobacillus", 
#             "Bacteroides", "Bacteroides", "Helicobacter"),
#   Species = c("Pseudomonas aeruginosa", "Pseudomonas aeruginosa", "Pseudomonas aeruginosa", 
#               "Escherichia coli", "Escherichia coli", "Escherichia coli", 
#               "Clostridium difficile", "Clostridium difficile", "Clostridium difficile", 
#               "Salmonella enterica", "Enterococcus faecalis", "Staphylococcus aureus", 
#               "Listeria monocytogenes", "Bacillus subtilis", "Lactobacillus plantarum", 
#               "Bacteroides fragilis", "Bacteroides vulgatus", "Helicobacter pylori")
# )
# 
# # Convert to matrix 
# taxa_matrix <- as.matrix(taxa_data[, -1])
# rownames(taxa_matrix) <- taxa_data$OTUID
# 
# # Step 2: Create an OTU table
# otu_data <- round(matrix(
#   c(5.1, 2.3, 1.5,    # Pseudomonas aeruginosa ASV1, ASV2, ASV3
#     12.4, 6.8, 5.9,   # Escherichia coli ASV4, ASV5, ASV6
#     15.2, 7.3, 6.9,   # Clostridium difficile ASV7, ASV8, ASV9
#     12.7, 19, 17.3,   # Salmonella enterica
#     21.4, 10.3, 14.6, # Enterococcus faecalis
#     13.1, 9.8, 4.6,   # Staphylococcus aureus
#     2.5, 1.8, 11.2,   # Listeria monocytogenes
#     7.1, 6.3, 12.7,   # Bacillus subtilis
#     5.8, 5.2, 11.4,   # Lactobacillus plantarum
#     17.1, 15.6, 19.2, # Bacteroides fragilis
#     12.7, 13.8, 12.5, # Bacteroides vulgatus
#     8.3, 7.8, 3.9),   # Helicobacter pylori
#   nrow = 18, ncol = 12, byrow = TRUE,   
#   dimnames = list(taxa_data$OTUID, paste0("Sample", 1:12))
# ))
# 
# # Step 3: Create sample metadata with 12 samples and 4 rep 
# sample_data <- data.frame(
#   SampleID = paste0("Sample", 1:12),  
#   Category = rep(c("Control", "Extreme Environment", "Normal Condition"), each = 4),  # 4 replicates for each condition
#   row.names = paste0("Sample", 1:12)  
# )
# 
# # Step 4: build the phyloseq 
# otu_table_ps <- otu_table(otu_data, taxa_are_rows = TRUE)
# tax_table_ps <- tax_table(taxa_matrix)
# sample_data_ps <- sample_data(sample_data)
# 
# physeq <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
# 
# # tidy up
# physeq<- tidy_phyloseq(physeq)
# 
# spiked_species_list <- list(
#   c("Pseudomonas aeruginosa"),
#   c("Escherichia coli"),
#   c("Clostridium difficile")
# )
# 
# spiked_cells_list <- c(10000, 20000, 15000)
# 
# # Step 6: Calculate the scaling factors after merging the redundant spikein species
# scaling_factors <- calculate_list_average_scaling_factors(merged_physeq_sum,
# spiked_species_list, spiked_cells_list, merge_method = "sum") # or max
# # Print the scaling factors for each OTU
# print(scaling_factors)
