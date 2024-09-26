#' Pre-process species in a phyloseq object by merging ASVs
#'
#' This function pre-processes species in a phyloseq object by merging ASVs based on a specified method (sum or max).
#' It retains the Genus and Species information for the merged taxa.
#'
#' @param physeq A `phyloseq::phyloseq` object containing the microbiome data.
#' @param species_name A character vector of species names to be processed.
#' @param merge_method The method to use for merging ASVs: `"sum"` or `"max"`. Default is `"sum"`.
#' @param output_file Optional. A file path to save the merged phyloseq object.
#' @return A `phyloseq::phyloseq` object with the specified species pre-processed.
#' @examples
#' \dontrun{
#' # Example for 16S rRNA data:
#' # ps is a phyloseq object
#' species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
#' 
#' # Merging species by summing or taking max ASV abundances
#' merged_physeq_sum <- Pre_processing_species(ps, species_name, 
#' merge_method = "sum", output_file = "merged_physeq_sum.rds")
#' 
#' # Example for ITS rDNA data:
#' # Subsetting Dekkera species
#' Dekkera <- phyloseq::subset_taxa(ps, Species == "Dekkera_bruxellensis")
#' species_name <- "Dekkera_bruxellensis"
#' 
#' # Merging species by summing ASV abundances
#' merged_physeq_sum_ITS <- Pre_processing_species(Dekkera, species_name, 
#' merge_method = "sum", output_file = "merged_physeq_sum_ITS.rds")
#' }
#' @importFrom phyloseq otu_table tax_table merge_taxa subset_taxa taxa_sums
#' @importFrom methods is
#' @export
Pre_processing_species <- function(physeq, species_name, merge_method = c("sum", "max"), output_file = NULL) {
  
  merge_method <- match.arg(merge_method)
  message("Starting pre-processing...")
  
  for (species in species_name) {
    message("Processing species: ", species)
    
    # Get the taxonomy table
    taxa_table <- as.data.frame(phyloseq::tax_table(physeq))
    
    # Find ASVs belonging to the specified species
    species_asvs <- rownames(taxa_table[taxa_table$Species == species, ])
    
    # Debugging: Print the ASVs found for the species
    message("ASVs found for species ", species, ": ", toString(species_asvs))
    
    if (length(species_asvs) > 1) {
      if (merge_method == "sum") {
        # Sum the abundances of all ASVs of the specified species
        sum_abundances <- colSums(phyloseq::otu_table(physeq)[species_asvs, , drop = FALSE])
        
        # Create a new OTU table with summed abundances
        new_otu_table <- phyloseq::otu_table(physeq)
        new_otu_table[species_asvs[1], ] <- sum_abundances
        new_otu_table <- new_otu_table[-which(rownames(new_otu_table) %in% species_asvs[-1]), ]
        
        # Ensure OTU table dimensions and names are correct
        new_otu_table <- phyloseq::otu_table(new_otu_table, taxa_are_rows = TRUE)
        
        # Update the phyloseq object with the new OTU table
        physeq@otu_table <- new_otu_table
        
        # Update the taxonomy table
        new_tax_table <- phyloseq::tax_table(physeq)
        new_tax_table <- new_tax_table[-which(rownames(new_tax_table) %in% species_asvs[-1]), ]
        
        # Ensure taxonomy table dimensions and names are correct
        new_tax_table <- phyloseq::tax_table(new_tax_table)
        
        # Retain the Genus and Species information
        phyloseq::tax_table(physeq)[species_asvs[1], "Genus"] <- taxa_table[species_asvs[1], "Genus"]
        phyloseq::tax_table(physeq)[species_asvs[1], "Species"] <- species
        
        message("Merged phyloseq object by summing abundances for species: ", species)
        
      } else if (merge_method == "max") {
        # Find the ASV with the highest abundance
        max_abundance_asv <- which.max(phyloseq::taxa_sums(physeq)[species_asvs])
        
        # Merge all ASVs of the specified species into the one with the highest abundance
        archetype <- species_asvs[max_abundance_asv]
        physeq <- phyloseq::merge_taxa(physeq, eqtaxa = species_asvs, archetype = archetype)
        
        # Retain the Genus and Species information
        phyloseq::tax_table(physeq)[archetype, "Genus"] <- taxa_table[archetype, "Genus"]
        phyloseq::tax_table(physeq)[archetype, "Species"] <- species
        
        message("Merged phyloseq object by selecting maximum abundances for species: ", species)
      }
      
      # Save the merged phyloseq object if output file provided
      if (!is.null(output_file)) {
        saveRDS(physeq, file = output_file)
        message("Merged phyloseq object saved to: ", output_file)
      }
    } else {
      message("No need to merge; only one ASV for species: ", species)
    }
  }
  
  # Ensure the dimensions and row names of OTU and taxonomy tables match
  if (!all(rownames(phyloseq::otu_table(physeq)) == rownames(phyloseq::tax_table(physeq)))) {
    stop("Mismatch between OTU and taxonomy table row names.")
  }
  
  # Additional debugging: print the dimensions and row names
  message("OTU table dimensions: ", toString(dim(phyloseq::otu_table(physeq))))
  message("Taxonomy table dimensions: ", toString(dim(phyloseq::tax_table(physeq))))
  message("OTU table row names: ", toString(head(rownames(phyloseq::otu_table(physeq)))))
  message("Taxonomy table row names: ", toString(head(rownames(phyloseq::tax_table(physeq)))))
  
  message("Pre-processing complete.")
  
  return(physeq)
}

# Example:
# 16S rRNA
# presence of 'spiked.volume' column in metadata
# spiked_cells <-1847
# species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
# merged_spiked_species<-"Tetragenococcus_halophilus"
# Tetragenococcus <- phyloseq::subset_taxa(physeq_16SASV,
# Species=="Tetragenococcus_halophilus" | Species=="Tetragenococcus_sp")
# hashcodes <- row.names(phyloseq::tax_table(Tetragenococcus))
# 
# # ITS rDNA
# # presence of 'spiked.volume' column in metadata
# spiked_cells <- 733
# species_name <- spiked_species<-merged_spiked_species<-"Dekkera_bruxellensis"
# Dekkera <- phyloseq::subset_taxa(physeq_ITSASV, Species=="Dekkera_bruxellensis")
# hashcodes <- row.names(phyloseq::tax_table(Dekkera))
# 
# 
# merged_physeq_sum <- Pre_processing_species(physeq_16SASV, species_name,
# merge_method = "sum", output_file = "merged_physeq_sum.rds")
# merged_physeq_max <- Pre_processing_species(physeq_16SASV, species_name,
# merge_method = "max", output_file = "merged_physeq_max.rds")
