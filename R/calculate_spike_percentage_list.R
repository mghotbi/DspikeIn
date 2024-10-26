#' Calculate Spike Percentage for List of Specified Taxa in a Phyloseq Object
#'
#' This function calculates the percentage of reads from specified spiked species in a \code{phyloseq} object.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as DOCX and CSV files.
#'
#' @param physeq A \code{phyloseq} object containing microbial data.
#' @param merged_spiked_species A character vector of spiked species to check in the phyloseq object.
#' @param output_path A character string specifying the path to save the output files. Default is \code{"merged_data.docx"}.
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is \code{c(0.1, 11)}.
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @examples
#' \dontrun{
#' # Example usage:
#' spiked_species_list <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#' result <- calculate_spike_percentage_list(physeq, merged_spiked_species = spiked_species_list, passed_range = c(0.1, 10))
#' print(result)
#' }
#' @importFrom phyloseq tax_table sample_names sample_sums otu_table merge_taxa prune_taxa
#' @importFrom dplyr left_join
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom utils write.csv
#' @export
calculate_spike_percentage_list <- function(physeq, merged_spiked_species = NULL, output_path = "merged_data.docx", passed_range = c(0.1, 11)) {
  
  # Validate the passed_range parameter
  if (!is.numeric(passed_range) || length(passed_range) != 2) {
    stop("passed_range must be a numeric vector of length 2.")
  }
  if (is.null(merged_spiked_species)) stop("You must provide 'merged_spiked_species'.")
  
  # Clean and flatten the species names
  cleaned_spiked_species <- trimws(unlist(merged_spiked_species))
  
  # Extract species names from tax_table and create a logical index for spiked species
  species_names <- phyloseq::tax_table(physeq)[, "Species"]
  spiked_taxa_idx <- species_names %in% cleaned_spiked_species
  
  # Subset the phyloseq object for only the spiked taxa
  spiked_taxa <- phyloseq::prune_taxa(spiked_taxa_idx, physeq)
  
  # Check if any spiked taxa are present
  if (phyloseq::ntaxa(spiked_taxa) == 0) stop("No samples contain the specified spiked taxa.")
  
  # Merge ASVs for the spiked taxa
  merged_spiked <- phyloseq::merge_taxa(spiked_taxa, phyloseq::taxa_names(spiked_taxa))
  
  # Align OTU and taxonomy tables
  merged_spiked <- fix_phyloseq_dimensions(merged_spiked)
  
  # Calculate total reads for each sample
  total_reads <- data.frame(Sample = phyloseq::sample_names(physeq), 
                            Total_Reads = phyloseq::sample_sums(phyloseq::otu_table(physeq)))
  
  # Calculate reads specific to the merged spiked taxa
  spiked_reads <- data.frame(Sample = phyloseq::sample_names(merged_spiked),
                             Total_Reads_spiked = phyloseq::sample_sums(phyloseq::otu_table(merged_spiked)))
  
  # Merge total reads and spiked reads into a single data frame using dplyr::left_join
  merged_data <- dplyr::left_join(total_reads, spiked_reads, by = "Sample")
  
  # Calculate the percentage of spiked taxa reads relative to total reads
  merged_data$Percentage <- (merged_data$Total_Reads_spiked / merged_data$Total_Reads) * 100
  
  # Categorize results as "passed" or "failed" based on the passed_range
  merged_data$Result <- ifelse(merged_data$Percentage >= passed_range[1] & merged_data$Percentage <= passed_range[2], "passed", "failed")
  
  # Create a formatted table using flextable
  ft <- flextable::flextable(merged_data) %>%
    flextable::fontsize(size = 10) %>%
    flextable::font(part = "all", fontname = "Inconsolata") %>%
    flextable::color(part = "header", color = "red4") %>%
    flextable::bold(part = "header") %>%
    flextable::italic()
  
  # Save the flextable as a Word document
  flextable::save_as_docx(ft, path = output_path)
  
  # Save merged data frame as CSV
  csv_path <- sub(".docx", ".csv", output_path)
  utils::write.csv(merged_data, file = csv_path, row.names = FALSE)
  
  # Print file save locations
  cat("Table saved in DOCX format:", output_path, "\n")
  cat("Merged data saved as CSV:", csv_path, "\n")
  
  # Return the merged data frame
  return(merged_data)
}

# Helper function to fix phyloseq table dimensions
#' @keywords internal
fix_phyloseq_dimensions <- function(physeq) {
  otu_table_names <- rownames(phyloseq::otu_table(physeq))
  tax_table_names <- rownames(phyloseq::tax_table(physeq))
  
  # Align OTU and taxonomy tables
  if (!all(otu_table_names %in% tax_table_names)) {
    physeq@tax_table <- phyloseq::tax_table(phyloseq::tax_table(physeq)[otu_table_names, , drop = FALSE])
  }
  
  return(physeq)
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
# spiked_species_list <- list(
#   c("Pseudomonas aeruginosa"),
#   c("Escherichia coli"),
#   c("Clostridium difficile")
# )
# 
# # Call the function to calculate the spike percentages
# result <- calculate_spike_percentage_list(merged_physeq_sum, merged_spiked_species = spiked_species_list, passed_range = c(0.1, 20))
# 
# # Print the result
# print(result)
