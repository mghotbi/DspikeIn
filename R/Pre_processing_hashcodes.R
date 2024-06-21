#' Pre-process phyloseq object based on hashcodes
#'
#' This function pre-processes a phyloseq object by subsetting, merging, and saving taxa
#' based on specified hashcodes and a merge method (sum or max). It retains the taxa 
#' information and creates various intermediate and final datasets for further processing.
#'
#' @param spiked_16S A phyloseq object containing the spiked dataset.
#' @param hashcodes A character vector of hashcodes (row names of the OTU table) to be processed.
#' @param merge_method The method to use for merging taxa: "sum" or "max".
#' @param output_prefix A character string to be used as a prefix for the output file names.
#' @return A phyloseq object with the processed data.
#' @examples
#' Tetra <- subset_taxa(physeqASV16, Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp")
#' hashcodes <- row.names(phyloseq::otu_table(Tetra))
#' processed_data_sum <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "sum", output_prefix = "merged_physeq_sum")
#' processed_data_max <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "max", output_prefix = "merged_physeq_max")
#' summ_count_phyloseq(processed_data_sum)
#' @export
Pre_processing_hashcodes <- function(spiked_16S, hashcodes, merge_method = c("sum", "max"), output_prefix = "merged_physeq") {
  merge_method <- match.arg(merge_method)
  message("Starting pre-processing...")
  
  # Check if input data and parameters are valid
  if (is.null(spiked_16S)) {
    stop("Error: spiked_16S dataset is NULL.")
  }
  if (!all(hashcodes %in% phyloseq::taxa_names(spiked_16S))) {
    stop("Error: One or more hashcodes not found in the dataset.")
  }
  
  if (merge_method == "sum") {
    # Use merge_taxa to merge taxa by summing their abundance
    spiked_16S_merged <- phyloseq::merge_taxa(spiked_16S, hashcodes)
  } else if (merge_method == "max") {
    message("Merging taxa using the 'max' method...")
    
    # Extract the OTU table and find the hashcode with the maximum abundance for each sample
    otu_tab <- phyloseq::otu_table(spiked_16S)
    max_abundances <- apply(otu_tab[hashcodes, , drop = FALSE], 2, max)
    max_hashcode_indices <- apply(otu_tab[hashcodes, , drop = FALSE], 2, which.max)
    max_hashcodes <- hashcodes[max_hashcode_indices]
    
    # Create a new OTU table with the maximum abundances
    new_otu_tab <- otu_tab[max_hashcodes[1], , drop = FALSE]
    new_otu_tab[] <- max_abundances
    
    # Retain the taxonomic information for the merged taxon
    new_tax_tab <- phyloseq::tax_table(spiked_16S)[max_hashcodes[1], , drop = FALSE]
    
    # Create a new phyloseq object with the updated OTU and taxonomy tables
    max_phyloseq <- phyloseq::phyloseq(phyloseq::otu_table(new_otu_tab, taxa_are_rows = TRUE), 
                                       phyloseq::tax_table(new_tax_tab),
                                       phyloseq::sample_data(spiked_16S))
    
    # Remove the hashcodes from the original phyloseq object
    spiked_16S_pruned <- phyloseq::prune_taxa(!phyloseq::taxa_names(spiked_16S) %in% hashcodes, spiked_16S)
    
    # Merge the pruned phyloseq object with the max abundance phyloseq object
    spiked_16S_merged <- phyloseq::merge_phyloseq(spiked_16S_pruned, max_phyloseq)
  }
  
  # Save the processed data
  processed_data_path <- paste0(output_prefix, "_processed.rds")
  saveRDS(spiked_16S_merged, processed_data_path)
  message("Saved processed data file: ", processed_data_path)
  
  message("Pre-processing complete.")
  return(spiked_16S_merged)
}

# Example usage:
# Tetra <- subset_taxa(physeqASV16, Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp")
# hashcodes <- row.names(phyloseq::otu_table(Tetra))
# processed_data_sum <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "sum", output_prefix = "merged_physeq_sum")
# processed_data_max <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "max", output_prefix = "merged_physeq_max")
