#' Calculate Spike Percentage for list of Specified Taxa in a Phyloseq Object
#'
#' This function calculates the percentage of reads from specified spiked species or hashcodes in a phyloseq object.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as a DOCX and CSV file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked species list to check in the phyloseq object. Default is NULL.
#' @param merged_spiked_hashcodes A character vector of spiked hashcodes list to check in the phyloseq object. Default is NULL.
#' @param output_path A character string specifying the path to save the output files. Default is "merged_data.docx".
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @examples
#' library(phyloseq)
#' # Load example data
#' physeq <- ps
#'
#' # Define the spiked species list
#' spiked_species_list <- list(
#'   c("Methylobacterium_phyllostachyos"),
#'   c("Methylorubrum_salsuginis"),
#'   c("Bosea_massiliensis"),
#'   c("Bacillus_decolorationis")
#' )
#'
#' # Use the spiked species list in the function
#' merged_spiked_species <- spiked_species_list
#' calculate_spike_percentage_list(ps, merged_spiked_species = merged_spiked_species, 
#' passed_range = c(0.1, 10))
#' @export
calculate_spike_percentage_list <- function(physeq, merged_spiked_species = NULL, merged_spiked_hashcodes = NULL, output_path = "merged_data.docx", passed_range = c(0.1, 11)) {
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("flextable", quietly = TRUE)) {
    stop("Package 'flextable' is required but not installed.")
  }
  
  library(phyloseq)
  library(dplyr)
  library(flextable)
  
  # Determine the taxonomic identifiers to use
  if (!is.null(merged_spiked_species)) {
    spiked_taxa <- subset_taxa(physeq, phyloseq::tax_table(physeq)[, "Species"] %in% unlist(merged_spiked_species))
  } else if (!is.null(merged_spiked_hashcodes)) {
    spiked_taxa <- subset_taxa(physeq, rownames(phyloseq::tax_table(physeq)) %in% merged_spiked_hashcodes)
  } else {
    stop("You must provide either 'merged_spiked_species' or 'merged_spiked_hashcodes'.")
  }
  
  # Check if there are any samples containing the spiked taxa
  if (ntaxa(spiked_taxa) == 0) {
    stop("No samples contain the specified spiked taxa.")
  }
  
  # Calculate total reads for samples
  total_reads <- data.frame(Sample = sample_names(physeq), 
                            Total_Reads = sample_sums(phyloseq::otu_table(physeq)))
  
  # Merge all ASVs rooted from spiked taxa into one ASV
  merged_spiked <- merge_taxa(spiked_taxa, taxa_names(spiked_taxa))
  
  # Calculate reads specific to merged spiked taxa
  spiked_reads <- data.frame(Sample = sample_names(merged_spiked),
                             Total_Reads = sample_sums(phyloseq::otu_table(merged_spiked)))
  
  # Merge total reads and spiked reads data frames
  merged_data <- merge(total_reads, spiked_reads, by = "Sample", suffixes = c("_total", "_spiked"))
  
  # Calculate the percentage of spiked taxa reads relative to total reads
  merged_data$Percentage <- (merged_data$Total_Reads_spiked / merged_data$Total_Reads_total) * 100
  
  # Categorize the results as "passed" or "failed" based on the passed_range
  merged_data$Result <- ifelse(merged_data$Percentage >= passed_range[1] & merged_data$Percentage <= passed_range[2], "passed", "failed")
  
  # Create flextable
  ft <- flextable(merged_data) %>% 
    flextable::fontsize(size = 10) %>% 
    flextable::font(part = "all", fontname = "Inconsolata") %>% 
    flextable::color(part = "header", color = "red4") %>% 
    flextable::bold(part = "header") %>% 
    flextable::italic() 
  
  # Save the flextable as a Word document
  save_as_docx(ft, path = output_path)
  
  # Save merged data frame as CSV
  csv_path <- sub(".docx", ".csv", output_path)
  write.csv(merged_data, file = csv_path, row.names = FALSE)
  
  # Print a message indicating the files have been saved
  cat("Table saved in docx format:", output_path, "\n")
  cat("Merged data saved as CSV:", csv_path, "\n")
  
  # Return the merged data frame
  return(merged_data)
}

# Example usage:
# Define the spiked species list
# spiked_species_list <- list(
#   c("Methylobacterium_phyllostachyos"),
#   c("Methylorubrum_salsuginis"),
#   c("Bosea_massiliensis"),
#   c("Bacillus_decolorationis")
# )
# #Use the spiked species list in the function
# merged_spiked_species <- spiked_species_list
# calculate_spike_percentage_list(ps, merged_spiked_species = merged_spiked_species, passed_range = c(0.1, 10))

# Define the spiked hashcodes
# merged_spiked_hashcodes <- c("hashcode1", "hashcode2")
# calculate_spike_percentage_list(physeq, merged_spiked_hashcodes = merged_spiked_hashcodes, passed_range = c(0.1, 10))

# Example of calculating summary stats
# result <- calculate_spike_percentage_list(spiked_ITS_OTU_scaled, merged_spiked_species = merged_spiked_species, passed_range = c(0.1, 35))
# calculate_summary_stats_table(result)
# result$Percentage
