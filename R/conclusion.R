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
    dplyr::summarize(
      mean_total_reads_spiked = mean(Total_Reads_spiked, na.rm = TRUE),
      sd_total_reads_spiked = sd(Total_Reads_spiked, na.rm = TRUE),
      se_total_reads_spiked = sd(Total_Reads_spiked, na.rm = TRUE) / sqrt(dplyr::n()),
      q25_total_reads_spiked = quantile(Total_Reads_spiked, 0.25, na.rm = TRUE),
      median_total_reads_spiked = median(Total_Reads_spiked, na.rm = TRUE),
      q75_total_reads_spiked = quantile(Total_Reads_spiked, 0.75, na.rm = TRUE),
      mean_percentage = mean(Percentage, na.rm = TRUE),
      sd_percentage = sd(Percentage, na.rm = TRUE),
      se_percentage = sd(Percentage, na.rm = TRUE) / sqrt(dplyr::n()),
      q25_percentage = quantile(Percentage, 0.25, na.rm = TRUE),
      median_percentage = median(Percentage, na.rm = TRUE),
      q75_percentage = quantile(Percentage, 0.75, na.rm = TRUE),
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
