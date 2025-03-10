#' @title Calculate Spike Percentage for Specified Taxa in a Phyloseq or TSE Object
#'
#' @description This function calculates the percentage of reads from specified spiked species or hashcodes in a microbiome dataset.
#' It merges the spiked taxa into one ASV, calculates the percentage of reads, categorizes the results as passed or failed,
#' and saves the results as a DOCX and CSV file.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` (TSE) object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked species to check in the dataset. Default is NULL.
#' @param merged_spiked_hashcodes A character vector of spiked hashcodes to check in the dataset. Default is NULL.
#' @param output_path A character string specifying the path to save the output files. Default is "merged_data.docx".
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#' @return A data frame containing the percentage of spiked taxa reads and the pass/fail results.
#' @importFrom phyloseq subset_taxa tax_table ntaxa sample_names sample_sums merge_taxa taxa_names otu_table
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom dplyr filter mutate pull summarise group_by ungroup desc all_of left_join
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#'   spiked_cells <- 1847
#'   species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- "Tetragenococcus_halophilus"
#'
#'   # Load phyloseq object from the DspikeIn package
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Pre-process species data by merging using the sum method
#'   merged_sum <- Pre_processing_species(physeq_16SOTU,
#'                                        species_name = species_name,
#'                                        merge_method = "sum")
#'
#'   # Calculate spike percentage within the specified range
#'   calculate_spike_percentage(merged_sum,
#'                              merged_spiked_species = merged_spiked_species,
#'                              passed_range = c(0.1, 20))
#' }
#' @export
calculate_spike_percentage <- function(obj, merged_spiked_species = NULL, merged_spiked_hashcodes = NULL,
                                       output_path = "merged_data.docx", passed_range = c(0.1, 11)) {
  suppressMessages({

    # Validate input type
    if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
      stop("\U0000274C	Input 'obj' must be a phyloseq or TreeSummarizedExperiment object.")
    }

    # Extract OTU table, taxonomy, and sample data using accessors
    otu_table <- get_otu_table(obj)
    tax_data <- get_tax_table(obj)
    sample_data <- get_sample_data(obj)

    # Identify spiked taxa based on species or hashcodes
    if (!is.null(merged_spiked_species)) {
      spiked_taxa <- rownames(tax_data)[tax_data[, "Species"] %in% merged_spiked_species]
    } else if (!is.null(merged_spiked_hashcodes)) {
      spiked_taxa <- rownames(tax_data)[rownames(tax_data) %in% merged_spiked_hashcodes]
    } else {
      stop("\U0000274C	You must provide either 'merged_spiked_species' or 'merged_spiked_hashcodes'.")
    }

    # Check if spiked taxa exist
    if (length(spiked_taxa) == 0) {
      stop("\U0000274C	No matching spiked taxa found in the dataset.")
    }

    # Get total reads per sample
    total_reads <- data.frame(
      Sample = colnames(otu_table),
      Total_Reads = colSums(otu_table)
    )

    # Subset OTU table to only include spiked taxa
    spiked_otu_table <- otu_table[spiked_taxa, , drop = FALSE]

    # Merge all spiked taxa into one entry
    spiked_reads <- data.frame(
      Sample = colnames(spiked_otu_table),
      Spiked_Reads = colSums(spiked_otu_table)
    )

    # Merge total and spiked reads into one table
    merged_data <- dplyr::left_join(total_reads, spiked_reads, by = "Sample")

    # Calculate the percentage of spiked taxa reads relative to total reads
    merged_data <- dplyr::mutate(merged_data, Percentage = (Spiked_Reads / Total_Reads) * 100)

    # Categorize samples as "passed" or "failed" based on the passed_range
    merged_data <- dplyr::mutate(merged_data, Result = ifelse(Percentage >= passed_range[1] & Percentage <= passed_range[2], "passed", "failed"))

    # Create flextable for the report
    ft <- flextable::flextable(merged_data) %>%
      flextable::fontsize(size = 10) %>%
      flextable::font(part = "all", fontname = "Inconsolata") %>%
      flextable::color(part = "header", color = "#6A0572") %>%
      flextable::bold(part = "header") %>%
      flextable::italic()

    # Save the flextable as a Word document
    flextable::save_as_docx(ft, path = output_path)

    # Save merged data as CSV
    csv_path <- sub(".docx", ".csv", output_path)
    utils::write.csv(merged_data, file = csv_path, row.names = FALSE)

    # Print message indicating saved files
    cat("\U0001F4C2 Table saved in docx format:", output_path, "\n")
    cat("\U0001F4C2 Merged data saved as CSV:", csv_path, "\n")
  })

  # Return the merged data frame
  return(merged_data)
}

# Example usage:
# Define the spiked species
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# calculate_spike_percentage(merged_TSE_sum, merged_spiked_species,
# passed_range = c(0.1, 20))


# merged_spiked_species<-"Dekkera_bruxellensis"
# result <- calculate_spike_percentage(merged_sum, merged_spiked_species,
# passed_range = c(0.1, 35))
# calculate_summary_stats_table(result)
# result$Percentage

# Define the spiked hashcodes
# merged_spiked_hashcodes <- c("hashcode1", "hashcode2")
# calculate_spike_percentage(physeq, merged_spiked_hashcodes = merged_spiked_hashcodes,
# passed_range = c(0.1, 10))


