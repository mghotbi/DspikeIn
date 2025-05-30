#' @title Compute Summary Statistics for Spiked Species
#'
#' @description
#' Computes per-sample spike-in summary statistics from a microbiome object (`phyloseq` or `TSE`),
#' generates a spike-in success report using `calculate_spike_percentage()`, and returns both the raw data and
#' a formatted summary table (`flextable`). The function also attempts to extract and retain the phylogenetic tree if present.
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
#'
#' @examples
#' ## -------------------------------
#' ## Example 1: Using phyloseq object
#' ## -------------------------------
#' library(DspikeIn)
#' data("physeq_16SOTU", package = "DspikeIn")
#'
#' # Merge spike-in species
#' species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#' merged_sum <- Pre_processing_species(
#'   obj = physeq_16SOTU,
#'   species_name = species_name,
#'   merge_method = "sum"
#' )
#'
#' # Compute summary statistics
#' output_doc <- file.path(tempdir(), "summary_phyloseq.docx")
#'
#' results_physeq <- conclusion(
#'   obj = merged_sum,
#'   merged_spiked_species = "Tetragenococcus_halophilus",
#'   max_passed_range = 20,
#'   output_path = output_doc
#' )
#' print(results_physeq$summary_stats)
#'
#' ## -----------------------------------------------
#' ## Example 2: Using TreeSummarizedExperiment object
#' ## -----------------------------------------------
#' tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#' output_doc_tse <- file.path(tempdir(), "summary_tse.docx")
#' results_tse <- conclusion(
#'   obj = tse_16SOTU,
#'   merged_spiked_species = "Tetragenococcus_halophilus",
#'   max_passed_range = 20,
#'   output_path = output_doc_tse
#' )
#' print(results_tse$summary_stats)
#'
#' @export
conclusion <- function(obj, merged_spiked_species, max_passed_range = 11, output_path = "merged_data.docx") {
  # Validate object type
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Unsupported object type: must be `phyloseq` or `TreeSummarizedExperiment`.")
  }

  # Extract components depending on the object type
  if (inherits(obj, "phyloseq")) {
    otu_table <- tryCatch(phyloseq::otu_table(obj), error = function(e) NULL)
    tax_data <- tryCatch(phyloseq::tax_table(obj), error = function(e) NULL)
    metadata <- tryCatch(phyloseq::sample_data(obj), error = function(e) NULL)
    phy_tree <- tryCatch(phyloseq::phy_tree(obj), error = function(e) NULL)
  } else {
    otu_table <- tryCatch(assay(obj), error = function(e) NULL)
    tax_data <- tryCatch(rowData(obj), error = function(e) NULL)
    metadata <- tryCatch(as.data.frame(SummarizedExperiment::colData(obj)), error = function(e) NULL)
    phy_tree <- NULL
  }

  if (!"spiked.volume" %in% colnames(metadata)) {
    stop("'spiked.volume' column not found in metadata.")
  }

  spike_success_report <- withCallingHandlers(
    calculate_spike_percentage(
      obj,
      merged_spiked_species,
      passed_range = c(0.1, max_passed_range)
    ),
    message = function(m) invokeRestart("muffleMessage")
  )


  if (!is.data.frame(spike_success_report)) {
    message("Extracting table from Word document...")
    spike_success_report <- tryCatch(
      read_spike_report(output_path),
      error = function(e) {
        stop(sprintf("Failed to extract table from Word document: %s", e$message))
      }
    )
  }

  if (!is.data.frame(spike_success_report)) {
    stop("Error: `calculate_spike_percentage()` did not return a valid dataframe.")
  }

  # Check for required columns
  required_cols <- c("Sample", "Total_Reads", "Spiked_Reads", "Percentage", "Result")
  missing_cols <- setdiff(required_cols, colnames(spike_success_report))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in `spike_success_report`: %s", paste(missing_cols, collapse = ", ")))
  }

  # Compute summary stats
  summary_stats_df <- spike_success_report |>
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

  # Format as flextable
  summary_stats_flex <- flextable::flextable(summary_stats_df) |>
    flextable::theme_vanilla() |>
    flextable::set_caption(" Summary Statistics of Spiked Species") |>
    flextable::bold(part = "header") |>
    flextable::bg(part = "header", bg = "#006400") |>
    flextable::color(part = "header", color = "white") |>
    flextable::align(align = "center", part = "all") |>
    flextable::fontsize(size = 11, part = "body") |>
    flextable::fontsize(size = 12, part = "header") |>
    flextable::width(j = seq_len(ncol(summary_stats_df)), width = 1.5) |>
    flextable::set_table_properties(width = 0.8, layout = "autofit")

  # Return results
  return(list(
    summary_stats = summary_stats_flex,
    full_report = spike_success_report,
    phy_tree = phy_tree
  ))
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
