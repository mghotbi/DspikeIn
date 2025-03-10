#' @title Calculate Spike Percentage for Specified Taxa
#' @description Computes the percentage of reads from specified spiked species in a
#' `phyloseq` or `TreeSummarizedExperiment` object. The function merges selected taxa
#' into a single ASV, calculates the percentage of reads per sample, classifies the results
#' as "passed" or "failed" based on a predefined range, and saves the results in DOCX
#' and CSV formats.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @param merged_spiked_species A character vector of spiked species names to check in the object.
#' @param output_path A character string specifying the file path for saving output (DOCX format).
#' Default is `"merged_data.docx"`.
#' @param passed_range A numeric vector of length 2 specifying the range of percentages for
#' categorization as "passed". Default is `c(0.1, 11)`.
#'
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail classification.
#' The results are also saved as a DOCX table and a CSV file.
#'
#' @details This function extracts relevant data from `phyloseq` or `TreeSummarizedExperiment`
#' objects using accessor functions, calculates the percentage of reads contributed by the specified
#' spiked taxa, and classifies each sample as "passed" or "failed" based on the provided threshold range.
#'
#' @importFrom phyloseq tax_table sample_names sample_sums otu_table
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr left_join
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom utils write.csv
#' @importFrom rlang .data
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load phyloseq object
#'   data("physeq", package = "DspikeIn")
#'
#'   # Define a list of spiked species
#'   spiked_species_list <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#'
#'   # Calculate spike percentage within the specified range
#'   result <- calculate_spike_percentage_list(
#'     obj = physeq,
#'     merged_spiked_species = spiked_species_list,
#'     passed_range = c(0.1, 10)
#'   )
#'
#'   # Print the result
#'   print(result)
#' }
#' }
#'
#' @export
calculate_spike_percentage_list <- function(obj, merged_spiked_species = NULL,
                                            output_path = "merged_data.docx",
                                            passed_range = c(0.1, 11)) {

  # Step 1: Validate Input Parameters
  if (!is.numeric(passed_range) || length(passed_range) != 2) {
    stop("passed_range must be a numeric vector of length 2.")
  }
  if (is.null(merged_spiked_species)) {
    stop("You must provide 'merged_spiked_species'.")
  }

  # Step 2: Extract Data Using Accessor Functions
  otu_data <- get_otu_table(obj) # Extract OTU table (phyloseq or TSE)
  tax_data <- get_tax_table(obj) # Extract taxonomy table (phyloseq or TSE)
  sample_data <- get_sample_data(obj) # Extract sample metadata (phyloseq or TSE)

  # Ensure that the species column exists in taxonomy data
  if (!"Species" %in% colnames(tax_data)) {
    stop("The taxonomy table does not contain a 'Species' column.")
  }

  # Step 3: Identify Spiked Taxa in the Taxonomy Table
  cleaned_spiked_species <- trimws(unlist(merged_spiked_species))

  # Get species names from taxonomy table
  species_names <- tax_data[, "Species", drop = FALSE]

  # Identify the taxa matching spiked species
  spiked_taxa_idx <- rownames(species_names)[species_names$Species %in% cleaned_spiked_species]

  # Check if any spiked taxa are found
  if (length(spiked_taxa_idx) == 0) {
    stop("\U0000274C No spiked taxa found in the dataset.")
  }

  # Step 4: Subset OTU Table for Spiked Taxa
  spiked_otu_data <- otu_data[spiked_taxa_idx, , drop = FALSE]

  # Merge Spiked Taxa into One ASV
  merged_spiked_reads <- colSums(spiked_otu_data, na.rm = TRUE)

  # Step 5: Calculate Total Reads for Each Sample
  total_reads <- colSums(otu_data, na.rm = TRUE)

  # Step 6: Create Data Frame with Calculations
  merged_data <- data.frame(
    Sample = colnames(otu_data),
    Total_Reads = total_reads,
    Total_Reads_spiked = merged_spiked_reads,
    Percentage = (merged_spiked_reads / total_reads) * 100
  )

  # Categorize Results Based on Passed Range
  merged_data$Result <- ifelse(
    merged_data$Percentage >= passed_range[1] & merged_data$Percentage <= passed_range[2],
    "passed",
    "failed"
  )

  # Step 7: Format and Save as DOCX Using flextable
  ft <- flextable::flextable(merged_data) %>%
    flextable::fontsize(size = 10) %>%
    flextable::font(part = "all", fontname = "Inconsolata") %>%
    flextable::color(part = "header", color = "red4") %>%
    flextable::bold(part = "header") %>%
    flextable::italic()

  # Save the table as DOCX and CSV
  flextable::save_as_docx(ft, path = output_path)
  csv_path <- sub(".docx", ".csv", output_path)
  utils::write.csv(merged_data, file = csv_path, row.names = FALSE)

  # Print save locations
  message("\U0001F4BE Table saved in DOCX format: ", output_path)
  message("\U0001F4BE Merged data saved as CSV: ", csv_path)

  # Return the final data frame
  return(merged_data)
}


# Example usage:
#
# # tidy up
# physeq<- tidy_phyloseq_tse(physeq)
# spiked_species_list <- list(
#   c("Pseudomonas aeruginosa"),
#   c("Escherichia coli"),
#   c("Clostridium difficile") )
#
# # Call the function to calculate the spike percentages
# result <- calculate_spike_percentage_list(merged_physeq_sum,
# merged_spiked_species = spiked_species_list, passed_range = c(0.1, 20))
#
# # Print the result
# print(result)
