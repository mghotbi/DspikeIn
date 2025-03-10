#' @title Compute Summary Statistics for Spiked Species
#'
#' @description Processes a phylogenetic dataset, extracts components,
#' calculates spiked species statistics, and generates summary statistics in `flextable` format.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param merged_spiked_species A character vector of spiked species names.
#' @param max_passed_range Numeric, maximum threshold for passing spike percentage.
#' @param output_path Character, file path for the `.docx` output from `calculate_spike_percentage()`.
#'
#' @return A list containing:
#'   \item{summary_stats}{A `flextable` summary of the spike statistics.}
#'   \item{full_report}{The full spiked species report as a `data.frame`.}
#'   \item{phy_tree}{The phylogenetic tree (if available).}
#'
#' @importFrom dplyr summarize mutate rename_with select filter
#' @importFrom tidyr pivot_wider
#' @importFrom SummarizedExperiment colData
#' @importFrom phyloseq otu_table tax_table sample_data phy_tree
#' @importFrom flextable flextable theme_vanilla set_caption autofit bg color align fontsize bold width
#' @importFrom officer read_docx docx_summary
#' @examples
#' \donttest{
#'   data("physeq_16SOTU", package="DspikeIn")
#'   species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- c("Tetragenococcus_halophilus")
#'   merged_sum <- Pre_processing_species(physeq_16SOTU, species_name, merge_method = "sum")
#'   max_passed_range <- 20
#'   results <- conclusion(merged_sum, merged_spiked_species, max_passed_range)
#'   results$summary_stats
#' }
#' @export
conclusion <- function(obj, merged_spiked_species, max_passed_range = 11, output_path = "merged_data.docx") {
  suppressMessages({

    #  Check if object is valid
    if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
      stop("\U0000274C Unsupported object type: must be `phyloseq` or `TreeSummarizedExperiment`.")
    }

    # Extract components depending on the object type
    if (inherits(obj, "phyloseq")) {
      otu_table <- tryCatch(phyloseq::otu_table(obj), error = function(e) NULL)
      tax_data <- tryCatch(phyloseq::tax_table(obj), error = function(e) NULL)
      metadata <- tryCatch(phyloseq::sample_data(obj), error = function(e) NULL)
      phy_tree <- tryCatch(phyloseq::phy_tree(obj), error = function(e) NULL)
    } else if (inherits(obj, "TreeSummarizedExperiment")) {
      otu_table <- tryCatch(assay(obj), error = function(e) NULL)
      tax_data <- tryCatch(rowData(obj), error = function(e) NULL)
      metadata <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(e) NULL)
      phy_tree <- NULL  # TSE does not support phylogenetic trees
    }

    #  Ensure required metadata columns exist
    if (!"spiked.volume" %in% colnames(metadata)) {
      stop("\U0000274C 'spiked.volume' column not found in metadata.")
    }

    #  Run `calculate_spike_percentage()`
    spike_success_report <- calculate_spike_percentage(obj, merged_spiked_species, passed_range = c(0.1, max_passed_range))

    #  Extract from `.docx` if `calculate_spike_percentage()` fails
    if (!is.data.frame(spike_success_report)) {
      message("\U0001F4DD Extracting table from Word document...")
      spike_success_report <- tryCatch({
        read_spike_report(output_path)
      }, error = function(e) {
        stop("\U0000274C Failed to extract table from Word document: ", e$message)
      })
    }

    #  Ensure output is a valid dataframe
    if (!is.data.frame(spike_success_report)) {
      stop("\U0000274C Error: `calculate_spike_percentage()` did not return a valid dataframe.")
    }

    # Check for required columns
    required_cols <- c("Sample", "Total_Reads", "Spiked_Reads", "Percentage", "Result")
    missing_cols <- setdiff(required_cols, colnames(spike_success_report))
    if (length(missing_cols) > 0) {
      stop("\U0000274C Missing columns in `spike_success_report`: ", paste(missing_cols, collapse = ", "))
    }

    #  Compute summary statistics
    summary_stats_df <- spike_success_report %>%
      dplyr::summarize(
        mean_total_reads_spiked = mean(Total_Reads, na.rm = TRUE),
        sd_total_reads_spiked = sd(Total_Reads, na.rm = TRUE),
        median_total_reads_spiked = median(Total_Reads, na.rm = TRUE),
        mean_percentage = mean(Percentage, na.rm = TRUE),
        sd_percentage = sd(Percentage, na.rm = TRUE),
        median_percentage = median(Percentage, na.rm = TRUE),
        passed_count = sum(Result == "passed"),
        failed_count = sum(Result == "failed")
      )

    # Convert summary statistics to `flextable`
    summary_stats_flex <- flextable::flextable(summary_stats_df) %>%
      flextable::theme_vanilla() %>%
      flextable::set_caption(" Summary Statistics of Spiked Species") %>%
      flextable::bold(part = "header") %>%
      flextable::bg(part = "header", bg = "#006400") %>%  # Dark Green Header
      flextable::color(part = "header", color = "white") %>%  # White text for contrast
      flextable::align(align = "center", part = "all") %>%
      flextable::fontsize(size = 11, part = "body") %>%
      flextable::fontsize(size = 12, part = "header") %>%
      flextable::width(j = 1:ncol(summary_stats_df), width = 1.5) %>%
      flextable::set_table_properties(width = 0.8, layout = "autofit")  # Prevents table overflow

    #  Return processed components
    return(list(
      summary_stats = summary_stats_flex,  # Flextable summary
      full_report = spike_success_report,  # Raw data as dataframe
      phy_tree = phy_tree  # Phylogenetic tree (if available)
    ))
  })
}


# Usage Example
# species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")

# merged_sum <- Pre_processing_species(physeq_16SOTU, species_name, merge_method = "sum")
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# max_passed_range <- 30

# results <- conclusion(merged_sum, merged_spiked_species, max_passed_range)
# print(results$summary_stats)
# head(results$full_report)

## Example usage:
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# max_passed_range <- 20
# output_path <- "spike_success_report.docx"

# Convert to absolute counts (example function)
# absolute <- convert_to_absolute_counts(physeq_16SOTU, scaling_factors)
# absolute$absolute_abundance_object
# physeq_absolute <- absolute$absolute_abundance_object

# # subset to exclude blanks
# physeq_adjusted <- phyloseq::subset_samples(physeq_absolute,
# sample.or.blank != "blank")
# absolute$summary_stats

## Run the conclusion function
# summary_stats <- conclusion(tse_16SOTU, merged_spiked_species, max_passed_range)

# Print the summary statistics
# print(summary_stats)
