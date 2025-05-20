#' @title Calculate Spike-in Percentage for Specified Taxa
#'
#' @description
#' Computes the percentage of reads attributed to specified spike-in taxa in a
#' `phyloseq` or `TreeSummarizedExperiment` object. The function merges spike-in taxa,
#' computes percentages, classifies samples into "passed" or "failed" based on a
#' user-defined threshold, and optionally exports DOCX and CSV reports.
#'
#' @details
#' The function automatically detects spike-in OTUs based on the `Species` column in
#' the taxonomy table. It works with both `phyloseq` and `TreeSummarizedExperiment`
#' objects and produces QC diagnostics commonly required for spike-in based absolute
#' quantification workflows.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param merged_spiked_species Character vector or list of spike-in species names (at the `Species` level).
#' @param output_path Optional. Character string specifying the output path (DOCX).
#' If `NULL` (default) results are not saved.
#' @param passed_range Numeric vector of length 2 specifying the accepted percentage range.
#' Default is `c(0.1, 11)`.
#'
#' @return A `data.frame` with:
#' \itemize{
#'   \item Sample
#'   \item Total_Reads
#'   \item Total_Reads_spiked
#'   \item Percentage (of spike-in reads)
#'   \item Result ("passed"/"failed")
#' }
#'
#' @section Notes:
#' - Assumes the taxonomy table contains a column named `Species`.
#' - Supports both `phyloseq` and `TreeSummarizedExperiment` objects.
#'
#' @importFrom phyloseq otu_table tax_table sample_data
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom utils write.csv
#' @importFrom stats na.omit
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq", package = "DspikeIn")
#'   spiked_species_list <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#'
#'   temp_docx <- file.path(tempdir(), "merged_result.docx")
#'   temp_csv <- sub(".docx", ".csv", temp_docx)
#'
#'   result <- calculate_spike_percentage_list(
#'     obj = physeq,
#'     merged_spiked_species = spiked_species_list,
#'     output_path = temp_docx,
#'     passed_range = c(0.1, 10)
#'   )
#'
#'   print(result)
#'
#'   # Clean up
#'   if (file.exists(temp_docx)) unlink(temp_docx, force = TRUE)
#'   if (file.exists(temp_csv)) unlink(temp_csv, force = TRUE)
#' }
#' @export
calculate_spike_percentage_list <- function(obj,
                                            merged_spiked_species,
                                            output_path = NULL,
                                            passed_range = c(0.1, 11)) {
  # --- Validate input ---
  if (!is.numeric(passed_range) || length(passed_range) != 2) {
    stop("'passed_range' must be a numeric vector of length 2.")
  }

  # Accept character() or list()
  merged_spiked_species <- unique(as.character(unlist(merged_spiked_species)))
  if (length(merged_spiked_species) == 0) {
    stop("You must provide at least one spike-in species.")
  }

  # --- Accessors ---
  otu_data <- if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::assay(obj)
  } else {
    as(phyloseq::otu_table(obj), "matrix")
  }

  tax_data <- if (inherits(obj, "TreeSummarizedExperiment")) {
    as.data.frame(SummarizedExperiment::rowData(obj))
  } else {
    as.data.frame(phyloseq::tax_table(obj))
  }

  # --- Check species column ---
  if (!"Species" %in% colnames(tax_data)) {
    stop("The taxonomy table must contain a 'Species' column.")
  }

  # --- Match spike-in taxa ---
  spikein_idx <- rownames(tax_data)[tax_data$Species %in% merged_spiked_species]

  if (length(spikein_idx) == 0) {
    stop("No matching spike-in species found in the dataset.")
  }

  # --- Calculate percentages ---
  spiked_otu_data <- otu_data[spikein_idx, , drop = FALSE]
  merged_spiked_reads <- colSums(spiked_otu_data, na.rm = TRUE)
  total_reads <- colSums(otu_data, na.rm = TRUE)

  merged_data <- data.frame(
    Sample = colnames(otu_data),
    Total_Reads = total_reads,
    Total_Reads_spiked = merged_spiked_reads,
    Percentage = (merged_spiked_reads / total_reads) * 100
  )

  merged_data$Result <- ifelse(
    merged_data$Percentage >= passed_range[1] & merged_data$Percentage <= passed_range[2],
    "passed", "failed"
  )

  # --- Save DOCX and CSV if output_path is provided ---
  if (!is.null(output_path)) {
    ft <- flextable::flextable(merged_data) |>
      flextable::fontsize(size = 10) |>
      flextable::font(part = "all", fontname = "Inconsolata") |>
      flextable::color(part = "header", color = "red4") |>
      flextable::bold(part = "header")

    flextable::save_as_docx(ft, path = output_path)
    csv_path <- sub("\\.docx$", ".csv", output_path)
    utils::write.csv(merged_data, file = csv_path, row.names = FALSE)

    message("\u2713 Results saved to: ", output_path)
    message("\u2713 CSV saved to: ", csv_path)
  }

  return(merged_data)
}


# Example usage:
#
# # tidy up
# physeq<- tidy_phyloseq_tse(physeq)
# spiked_species_list<- merged_spiked_species <- c("Pseudomonas aeruginosa",
# "Escherichia coli", "Clostridium difficile")

#
# # Call the function to calculate the spike percentages
# result <- calculate_spike_percentage_list(merged_physeq_sum,
# merged_spiked_species = spiked_species_list, passed_range = c(0.1, 20))
#
# # Print the result
# print(result)
