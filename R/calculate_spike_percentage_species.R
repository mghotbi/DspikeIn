#' Calculate Spike Percentage for Specified Species in Phyloseq Object
#'
#' This function calculates the percentage of reads from specified spiked species in a phyloseq object.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as a DOCX and CSV file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param spiked_species A character vector of spiked species to check in the phyloseq object.
#' @param identifier_type A character string specifying the type of identifier to use ("species" or "hashcode"). Default is "species".
#' @param output_path A character string specifying the path to save the output files. Default is NULL, which saves the report as "merged_data.docx".
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @examples
#' # Example usage:
#' # Define the spiked species
#' spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'
#' # Calculate spike percentage and generate the report
#' calculate_spike_percentage_species(spiked_16S, spiked_species, identifier_type = "species", output_path = NULL, passed_range = c(0.1, 11))
#'
#' # Save the results
#' merged_data <- calculate_spike_percentage_species(spiked_16S, spiked_species)
#' @export
calculate_spike_percentage_species <- function(physeq, spiked_species, identifier_type = "species", output_path = NULL, passed_range = c(0.1, 11)) {
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
  if (identifier_type == "species") {
    spiked_taxa <- subset_taxa(physeq, Species %in% spiked_species)
  } else if (identifier_type == "hashcode") {
    spiked_taxa <- subset_taxa(physeq, row.names(tax_table(physeq)) %in% spiked_species)
  } else {
    stop("Invalid identifier_type. Use 'species' or 'hashcode'.")
  }
  
  # Check if there are any samples containing the spiked taxa
  if (ntaxa(spiked_taxa) == 0) {
    stop("No samples contain the spiked species.")
  }
  
  # Calculate total reads for 16S samples
  total_reads_16S <- data.frame(Sample = sample_names(physeq), 
                                Total_Reads = sample_sums(otu_table(physeq)))
  
  # Merge all ASVs rooted from spiked taxa into one ASV
  merged_spiked <- merge_taxa(spiked_taxa, taxa_names(spiked_taxa))
  
  # Calculate reads specific to merged spiked species
  spiked_reads <- data.frame(Sample = sample_names(merged_spiked),
                             Total_Reads = sample_sums(otu_table(merged_spiked)))
  
  # Merge total reads and spiked reads data frames
  merged_data <- merge(total_reads_16S, spiked_reads, by = "Sample", suffixes = c("_total", "_spiked"))
  
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
  
  # Set default output directory if none provided
  if (is.null(output_path)) {
    output_path <- "merged_data.docx"
  }
  
  # Save the flextable as a Word document
  save_as_docx(ft, path = output_path)
  
  # Save merged data frame as CSV
  csv_path <- sub(".docx", ".csv", output_path)
  write.csv(merged_data, file = csv_path, row.names = FALSE)
  
  # Print a message 
  cat("Table saved in docx format:", output_path, "\n")
  cat("Merged data saved as CSV:", csv_path, "\n")
  
  # Return the merged data frame
  return(merged_data)
}

# Example:
# Define the spiked species
#spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#spiked_species <-"Dekkera_bruxellensis"
# Calculate spike percentage and generate the report
#calculate_spike_percentage_species(spiked_16S, spiked_species, identifier_type = "species", output_path = NULL, passed_range = c(0.1, 11))
#merged_data <- calculate_spike_percentage_species(spiked_16S, spiked_species)
