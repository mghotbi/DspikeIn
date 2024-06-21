#' Calculate Spike Percentage for Specified Taxa in a Phyloseq Object
#'
#' This function calculates the percentage of reads from specified spiked species or hashcodes in a phyloseq object.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as a DOCX and CSV file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked species to check in the phyloseq object. Default is NULL.
#' @param merged_spiked_hashcodes A character vector of spiked hashcodes to check in the phyloseq object. Default is NULL.
#' @param output_path A character string specifying the path to save the output files. Default is NULL, which saves the report as "merged_data.docx".
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @examples
#' # Example usage:
#' # Define the spiked species
#' merged_spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#' calculate_spike_percentage(physeq, merged_spiked_species = merged_spiked_species)
#'
#' # Define the spiked hashcodes
#' merged_spiked_hashcodes <- c("hashcode1", "hashcode2")
#' calculate_spike_percentage(physeq, merged_spiked_hashcodes = merged_spiked_hashcodes)
#' @export
calculate_spike_percentage <- function(physeq, merged_spiked_species = NULL, merged_spiked_hashcodes = NULL, output_path = NULL, passed_range = c(0.1, 11)) {
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
    spiked_taxa <- subset_taxa(physeq, phyloseq::tax_table(physeq)[, "Species"] %in% merged_spiked_species)
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
  
  # Set default output directory if none provided
  if (is.null(output_path)) {
    output_path <- "merged_data.docx"
  }
  
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

#' Generate Spike Success Report and Summary Statistics
#'
#' This function calculates the spike success report for a given phyloseq object, 
#' summarizes the results, and saves the summary statistics in both DOCX and CSV formats.
#' It uses the `flextable` package to create and format the tables in the DOCX file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked species to check in the phyloseq object.
#' @param max_passed_range A numeric value specifying the maximum acceptable spike percentage. Default is 11.
#' @param output_path A character string specifying the path to save the output files. Default is NULL, which saves the report as "spike_success_report.docx".
#' @return A data frame containing the summary statistics of the spike success report.
#' @examples
#' # Define the parameters
#' merged_spiked_species <- c("Dekkera_bruxellensis")
#' max_passed_range <- 11
#' output_path <- "spike_success_report.docx"
#'
#' # Subset the phyloseq object to exclude blanks
#' physeq_ITS_adj_scaled_perc <- subset_samples(physeq_ITS_OTU, sample_or_blank != "blank")
#'
#' # Generate the spike success report and summary statistics
#' conclusion(physeq_ITS_adj_scaled_perc, merged_spiked_species, max_passed_range, output_path)
#' @export
conclusion <- function(physeq, merged_spiked_species, max_passed_range = 11, output_path = NULL) {
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
  
  # Calculate spike percentage and generate spike success report
  spike_success_report <- calculate_spike_percentage(physeq, merged_spiked_species = merged_spiked_species, passed_range = c(0.1, max_passed_range), output_path = output_path)
  
  print("Structure of spike_success_report:")
  print(str(spike_success_report))
  
  # Debug: Print the first few rows of spike_success_report
  print("First few rows of spike_success_report:")
  print(head(spike_success_report))
  
  # Check if 'Result' column exists
  if (!"Result" %in% colnames(spike_success_report)) {
    stop("The 'Result' column is not present in the spike success report. Please check the calculate_spike_percentage function.")
  }
  
  # Convert the 'Result' column to lowercase
  spike_success_report$Result <- tolower(spike_success_report$Result)
  
  # Remove rows with NA in the 'Result' column
  spike_success_report <- spike_success_report %>% dplyr::filter(!is.na(Result))
  
  # Debug: Print unique values in the 'Result' column after filtering and converting to lowercase
  print("Unique values in the 'Result' column after filtering and converting to lowercase:")
  print(unique(spike_success_report$Result))
  
  # Calculate summary statistics
  summary_stats <- spike_success_report %>%
    summarise(across(where(is.numeric), list(
      mean = ~mean(., na.rm = TRUE),
      sd = ~ifelse(all(is.na(.)), NA, sd(., na.rm = TRUE)),
      se = ~ifelse(all(is.na(.)), NA, sd(., na.rm = TRUE) / sqrt(length(.))),
      q25 = ~quantile(., 0.25, na.rm = TRUE),
      median = ~median(., na.rm = TRUE),
      q75 = ~quantile(., 0.75, na.rm = TRUE)
    )),
    passed_count = sum(Result == "passed"),
    failed_count = sum(Result == "failed")
    )
  
  # Create flextable
  ft <- flextable::flextable(summary_stats) %>% 
    flextable::fontsize(size = 10) %>% 
    flextable::font(part = "all", fontname = "Inconsolata") %>% 
    flextable::color(part = "header", color = "red4") %>% 
    flextable::bold(part = "header") %>% 
    flextable::italic() 
  
  # Set default output directory if none provided
  if (is.null(output_path)) {
    output_path <- "spike_success_report.docx"
  }
  
  # Save the flextable as a Word document
  flextable::save_as_docx(ft, path = output_path)
  
  # Save summary statistics as CSV
  csv_path <- sub(".docx", ".csv", output_path)
  write.csv(summary_stats, file = csv_path, row.names = FALSE)
  
  # Print messages
  cat("Table saved in docx format:", output_path, "\n")
  cat("Summary statistics saved as CSV:", csv_path, "\n")
  
  # Return the summary statistics
  return(summary_stats)
}

# Example usage:
# Define the parameters
#merged_spiked_species <- c("Tetragenococcus_halophilus")
#max_passed_range <- 35

# Subset the phyloseq object to exclude blanks/optional
#physeq_16S_adj_scaled_perc <- subset_samples(Spiked_16S_ASV_scaled, sample.or.blank != "blank")

# Generate the spike success report and summary statistics
#summary_stats <- conclusion(physeq_16S_adj_scaled_perc, merged_spiked_species, max_passed_range)
#print(summary_stats)
