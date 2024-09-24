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
#' \dontrun{
#' # Define the parameters
#' merged_spiked_species <- c("Tetragenococcus_halophilus") 
#' max_passed_range <- 12 
#' output_path <- "spike_success_report.docx"
#' 
#' # Convert the phyloseq object to absolute counts using scaling factors
#' absolute <- convert_to_absolute_counts(merged_physeq_sum, scaling_factors) 
#' absolute_counts <- absolute$absolute_counts 
#' physeq_absolute <- absolute$physeq_obj 
#' 
#' # Subset the phyloseq object to exclude blank samples (optional step)
#' physeq_16S_adj_scaled_perc <- phyloseq::subset_samples(physeq_absolute, sample.or.blank != "blank")
#' 
#' # Generate the spike success report and calculate summary statistics
#' summary_stats <- conclusion(physeq_16S_adj_scaled_perc, 
#'                             merged_spiked_species, 
#'                             max_passed_range, 
#'                             output_path)
#' 
#' # Print the summary statistics
#' print(summary_stats)
#' }
#' @importFrom phyloseq subset_samples
#' @importFrom dplyr filter summarize n
#' @importFrom stats sd quantile median
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom utils write.csv
#' @export
conclusion <- function(physeq, merged_spiked_species, max_passed_range = 11, output_path = NULL) {
  suppressMessages({
    # Calculate spike percentage and generate spike success report
    spike_success_report <- calculate_spike_percentage(physeq, 
                                                       merged_spiked_species = merged_spiked_species, 
                                                       passed_range = c(0.1, max_passed_range), 
                                                       output_path = output_path)
    
    # Print structure of spike_success_report for debugging
    cat("Structure of spike_success_report:\n")
    str(spike_success_report)
    
    # Debug: Print the first few rows of spike_success_report
    cat("First few rows of spike_success_report:\n")
    print(head(spike_success_report))
    
    # Check if 'Result' column exists
    if (!"Result" %in% colnames(spike_success_report)) {
      stop("The 'Result' column is not present in the spike success report. Please check the calculate_spike_percentage function.")
    }
    
    # Convert the 'Result' column to lowercase
    spike_success_report$Result <- tolower(spike_success_report$Result)
    
    # Remove rows with NA in the 'Result' column
    spike_success_report <- dplyr::filter(spike_success_report, !is.na(Result))
    
    # Debug: Print unique values in the 'Result' column 
    cat("Unique values in the 'Result' column after filtering :\n")
    print(unique(spike_success_report$Result))
    
    # Calculate summary statistics
    summary_stats <- dplyr::summarize(
      spike_success_report,
      mean_total_reads_spiked = mean(Total_Reads_spiked, na.rm = TRUE),
      sd_total_reads_spiked = stats::sd(Total_Reads_spiked, na.rm = TRUE),
      se_total_reads_spiked = stats::sd(Total_Reads_spiked, na.rm = TRUE) / sqrt(dplyr::n()),
      q25_total_reads_spiked = stats::quantile(Total_Reads_spiked, 0.25, na.rm = TRUE),
      median_total_reads_spiked = stats::median(Total_Reads_spiked, na.rm = TRUE),
      q75_total_reads_spiked = stats::quantile(Total_Reads_spiked, 0.75, na.rm = TRUE),
      mean_percentage = mean(Percentage, na.rm = TRUE),
      sd_percentage = stats::sd(Percentage, na.rm = TRUE),
      se_percentage = stats::sd(Percentage, na.rm = TRUE) / sqrt(dplyr::n()),
      q25_percentage = stats::quantile(Percentage, 0.25, na.rm = TRUE),
      median_percentage = stats::median(Percentage, na.rm = TRUE),
      q75_percentage = stats::quantile(Percentage, 0.75, na.rm = TRUE),
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
    utils::write.csv(summary_stats, file = csv_path, row.names = FALSE)
    
    # Print messages
    cat("Table saved in docx format:", output_path, "\n")
    cat("Summary statistics saved as CSV:", csv_path, "\n")
    
    # Return the summary statistics
    return(summary_stats)
  })
}

# Example usage:
# Define the parameters
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# max_passed_range <- 35
# absolute <- convert_to_absolute_counts(merged_physeq_sum, scaling_factors)
# absolute_counts <- absolute$absolute_counts
# physeq_absolute <- absolute$physeq_obj
# # # Subset the phyloseq object to exclude blanks/optional
# physeq_16S_adj_scaled_perc <- phyloseq::subset_samples(physeq_absolute, sample.or.blank != "blank")
# #
# # Generate the spike success report and summary statistics
# summary_stats <- conclusion(physeq_16S_adj_scaled_perc, merged_spiked_species, max_passed_range)
# print(summary_stats)
