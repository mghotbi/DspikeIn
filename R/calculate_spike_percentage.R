#' Calculate Spike Percentage for Specified Taxa in a Phyloseq Object
#'
#' This function calculates the percentage of reads from specified spiked species or hashcodes in a phyloseq object.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as a DOCX and CSV file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked species to check in the phyloseq object. Default is NULL.
#' @param merged_spiked_hashcodes A character vector of spiked hashcodes to check in the phyloseq object. Default is NULL.
#' @param output_path A character string specifying the path to save the output files. Default is "merged_data.docx".
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @importFrom phyloseq subset_taxa tax_table ntaxa sample_names sample_sums merge_taxa taxa_names otu_table
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom dplyr filter mutate pull summarise group_by ungroup desc all_of left_join
#' @importFrom magrittr %>%
#' @export
#' @examples
#' \dontrun{
#' # Define the spiked species
#' merged_spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#' 
#' # Calculate spike percentage for the defined spiked species within the range 0.1 to 10
#' calculate_spike_percentage(physeq, 
#'                            merged_spiked_species = merged_spiked_species, 
#'                            passed_range = c(0.1, 10))
#' }

calculate_spike_percentage <- function(physeq, merged_spiked_species = NULL, merged_spiked_hashcodes = NULL, output_path = "merged_data.docx", passed_range = c(0.1, 11)) {
  suppressMessages({
    # Determine the taxonomic identifiers to use
    if (!is.null(merged_spiked_species)) {
      spiked_taxa <- phyloseq::subset_taxa(physeq, phyloseq::tax_table(physeq)[, "Species"] %in% merged_spiked_species)
    } else if (!is.null(merged_spiked_hashcodes)) {
      spiked_taxa <- phyloseq::subset_taxa(physeq, rownames(phyloseq::tax_table(physeq)) %in% merged_spiked_hashcodes)
    } else {
      stop("You must provide either 'merged_spiked_species' or 'merged_spiked_hashcodes'.")
    }
    
    # Check if there are any samples containing the spiked taxa
    if (phyloseq::ntaxa(spiked_taxa) == 0) {
      stop("No samples contain the specified spiked taxa.")
    }
    
    # Calculate total reads for samples
    total_reads <- data.frame(Sample = phyloseq::sample_names(physeq), 
                              Total_Reads = phyloseq::sample_sums(phyloseq::otu_table(physeq)))
    
    # Merge all ASVs rooted from spiked taxa into one ASV
    merged_spiked <- phyloseq::merge_taxa(spiked_taxa, phyloseq::taxa_names(spiked_taxa))
    
    # Calculate reads specific to merged spiked taxa
    spiked_reads <- data.frame(Sample = phyloseq::sample_names(merged_spiked),
                               Total_Reads = phyloseq::sample_sums(phyloseq::otu_table(merged_spiked)))
    
    # Merge total reads and spiked reads data frames (use suffixes if column names overlap)
    if ("Total_Reads" %in% colnames(total_reads) & "Total_Reads" %in% colnames(spiked_reads)) {
      merged_data <- dplyr::left_join(total_reads, spiked_reads, by = "Sample", suffix = c("_total", "_spiked"))
    } else {
      merged_data <- dplyr::left_join(total_reads, spiked_reads, by = "Sample")
    }
    
    # Calculate the percentage of spiked taxa reads relative to total reads
    merged_data <- dplyr::mutate(merged_data, Percentage = (Total_Reads_spiked / Total_Reads_total) * 100)
    
    # Categorize the results as "passed" or "failed" based on the passed_range
    merged_data <- dplyr::mutate(merged_data, Result = ifelse(Percentage >= passed_range[1] & Percentage <= passed_range[2], "passed", "failed"))
    
    # Create flextable
    ft <- flextable::flextable(merged_data) %>% 
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
    flextable::save_as_docx(ft, path = output_path)
    
    # Save merged data frame as CSV
    csv_path <- sub(".docx", ".csv", output_path)
    write.csv(merged_data, file = csv_path, row.names = FALSE)
    
    # Print a message indicating the files have been saved
    cat("Table saved in docx format:", output_path, "\n")
    cat("Merged data saved as CSV:", csv_path, "\n")
  })
  
  # Return the merged data frame
  return(merged_data)
}
# Example usage:
# Define the spiked species
# merged_spiked_species <- c("Tetragenococcus_halophilus","Tetragenococcus_sp")
# calculate_spike_percentage(physeq_ITSOTU, merged_spiked_species, 
# passed_range = c(0.1, 10))


# merged_spiked_species<-"Dekkera_bruxellensis"
# result <- calculate_spike_percentage(spiked_ITS_OTU_scaled, merged_spiked_species,
# passed_range = c(0.1, 35))
# calculate_summary_stats_table(result)
# result$Percentage

# Define the spiked hashcodes
# merged_spiked_hashcodes <- c("hashcode1", "hashcode2")
# calculate_spike_percentage(physeq, merged_spiked_hashcodes = merged_spiked_hashcodes,
# passed_range = c(0.1, 10))


