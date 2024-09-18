#' Pre-process a list of spiked-in species in a phyloseq object
#'
#' This function pre-processes a list of spiked-in species in a phyloseq object
#' by merging ASVs based on a specified method (sum or max). It retains the Genus
#' and Species information for the merged taxa as a single entry in the taxonomy table.
#'
#' @param physeq A \code{phyloseq} object containing the microbiome data.
#' @param spiked_species A character vector of spiked-in species names to be processed. Species names should match the \code{Genus Species} format in the phyloseq object (e.g., "Pseudomonas aeruginosa").
#' @param merge_method A character string specifying the method to use for merging ASVs. Either \code{"sum"} to sum the counts, or \code{"max"} to keep the maximum abundance ASV. Default is \code{"sum"}.
#' @param output_file Optional. A file path to save the merged phyloseq object. If provided, the merged object will be saved as an \code{.rds} file.
#' @return A \code{phyloseq} object with the spiked-in species pre-processed.
#' @examples
#' \dontrun{
#' # Example usage
#' spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#' merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
#' merged_physeq_sum@tax_table
#' merged_physeq_sum@otu_table
#' }
#' @importFrom phyloseq otu_table tax_table 
#' @export
Pre_processing_species_list <- function(physeq, spiked_species, merge_method = c("sum", "max"), output_file = NULL) {
  # Match the merge method argument
  merge_method <- match.arg(merge_method)
  
  # Suppress messages for loading required packages
  suppressMessages({
    requireNamespace("phyloseq", quietly = TRUE)
  })
  
  message("Starting pre-processing for spiked-in species...")
  
  # Get the taxonomy and OTU tables from the phyloseq object
  taxa_table <- as.data.frame(phyloseq::tax_table(physeq))
  otu_table_data <- as(phyloseq::otu_table(physeq), "matrix")
  
  # Loop over each species in the spiked_species list
  for (species in spiked_species) {
    message("Processing species: ", species)
    
    # Find all ASVs belonging to the current species (Genus + Species format)
    species_asvs <- rownames(taxa_table)[taxa_table$Species == species]
    
    if (length(species_asvs) > 1) {
      message("Merging ", length(species_asvs), " ASVs for species: ", species)
      
      # Merge the OTU counts based on the specified method
      if (merge_method == "sum") {
        merged_abundances <- colSums(otu_table_data[species_asvs, , drop = FALSE])
      } else if (merge_method == "max") {
        max_abundance_asv <- species_asvs[which.max(rowSums(otu_table_data[species_asvs, , drop = FALSE]))]
        merged_abundances <- otu_table_data[max_abundance_asv, ]
      }
      
      # Update the OTU table: retain the first ASV, update its counts, and remove others
      otu_table_data[species_asvs[1], ] <- merged_abundances
      otu_table_data <- otu_table_data[-which(rownames(otu_table_data) %in% species_asvs[-1]), ]
      
      # Update the taxonomy table: retain the first ASV's taxonomy and remove others
      taxa_table <- taxa_table[-which(rownames(taxa_table) %in% species_asvs[-1]), ]
      
      message("Successfully merged ASVs for species: ", species)
      
    } else if (length(species_asvs) == 1) {
      message("Only one ASV found for species: ", species, "; no merging required.")
    } else {
      message("No ASVs found for species: ", species)
    }
  }
  
  # Ensure that OTU and taxonomy tables have matching row names
  if (!all(rownames(otu_table_data) == rownames(taxa_table))) {
    stop("Mismatch between OTU and taxonomy table row names.")
  }
  
  # Update the phyloseq object with the modified OTU and taxonomy tables
  physeq@otu_table <- phyloseq::otu_table(as.matrix(otu_table_data), taxa_are_rows = TRUE)
  physeq@tax_table <- phyloseq::tax_table(as.matrix(taxa_table))
  
  # Optional: Save the merged phyloseq object if an output file path is provided
  if (!is.null(output_file)) {
    saveRDS(physeq, file = output_file)
    message("Merged phyloseq object saved to: ", output_file)
  }
  
  message("Pre-processing complete.")
  
  # Return the modified phyloseq object
  return(physeq)
}

# Example usage:
# Step 1: Create a taxonomy table
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
# # Example usage:
# spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
# merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
# merged_physeq_sum@tax_table
# merged_physeq_sum@otu_table
# 
# physeq <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
# 
# # tidy up
# physeq<- tidy_phyloseq(physeq)
# 
# spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
# merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
# merged_physeq_sum@tax_table
# merged_physeq_sum@otu_table
