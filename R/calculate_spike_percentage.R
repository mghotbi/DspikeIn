#' @title Calculate Spike Percentage for Specified Taxa in a Phyloseq or TSE Object
#'
#' @description This function calculates the percentage of reads from specified spiked species or hashcodes in a microbiome dataset.
#' It merges the spiked taxa into one ASV/OTU, calculates the percentage of reads, categorizes the results as passed or failed,
#' and optionally saves the results as DOCX and CSV files.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing the microbial data.
#' @param merged_spiked_species A character vector of spiked taxa names (can be from any taxonomic level).
#' @param merged_spiked_hashcodes A character vector of spiked hashcodes (ASV/OTU IDs) to check in the dataset. Default is NULL.
#' @param output_file A character string specifying the path to save the output files. Default is NULL (no files are written).
#' @param passed_range A numeric vector of length 2 specifying the range of percentages to categorize results as "passed". Default is c(0.1, 11).
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Sample}{Sample identifier.}
#'   \item{Total_Reads}{Total number of reads in the sample.}
#'   \item{Spiked_Reads}{Number of reads mapped to the spike-in taxa.}
#'   \item{Percentage}{Percentage of spike-in reads (Spiked_Reads / Total_Reads * 100).}
#'   \item{Result}{Quality control result, either \code{"passed"} or \code{"failed"}, based on the specified range.}
#' }
#'
#' @importFrom phyloseq otu_table tax_table
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom dplyr left_join mutate
#' @importFrom utils write.csv
#' @seealso \code{\link{Pre_processing_species}}, \code{\link{calculate_spike_percentage}}
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load example phyloseq object
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # ----------- Phyloseq Example -----------
#'   species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- "Tetragenococcus_halophilus"
#'
#'   # Pre-process the phyloseq object to merge spike-in taxa
#'   merged_physeq <- Pre_processing_species(
#'     physeq_16SOTU,
#'     species_name = species_name,
#'     merge_method = "sum"
#'   )
#'
#'   # Perform spike-in percentage calculation
#'   output_docx <- file.path(tempdir(), "spike_summary_physeq.docx")
#'   result_physeq <- calculate_spike_percentage(
#'     obj = merged_physeq,
#'     merged_spiked_species = merged_spiked_species,
#'     output_file = output_docx,
#'     passed_range = c(0.1, 20)
#'   )
#'   print(result_physeq)
#'
#'   # ----------- TreeSummarizedExperiment Example -----------
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   merged_tse <- Pre_processing_species(
#'     tse_16SOTU,
#'     species_name = species_name,
#'     merge_method = "sum"
#'   )
#'
#'   output_docx_tse <- file.path(tempdir(), "spike_summary_tse.docx")
#'   result_tse <- calculate_spike_percentage(
#'     obj = merged_tse,
#'     merged_spiked_species = merged_spiked_species,
#'     output_file = output_docx_tse,
#'     passed_range = c(0.1, 20)
#'   )
#'   print(result_tse)
#'
#'   # Clean up temporary files
#'   if (file.exists(output_docx)) unlink(output_docx, force = TRUE)
#'   if (file.exists(output_docx_tse)) unlink(output_docx_tse, force = TRUE)
#' }
#' @export
calculate_spike_percentage <- function(obj,
                                       merged_spiked_species = NULL,
                                       merged_spiked_hashcodes = NULL,
                                       output_file = NULL,
                                       passed_range = c(0.1, 11)) {
  # Validate input
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Input 'obj' must be a phyloseq or TreeSummarizedExperiment object.")
  }

  if (is.null(merged_spiked_species) && is.null(merged_spiked_hashcodes)) {
    stop("You must provide either 'merged_spiked_species' or 'merged_spiked_hashcodes'.")
  }

  if (!is.numeric(passed_range) || length(passed_range) != 2) {
    stop("'passed_range' must be a numeric vector of length 2.")
  }

  # Extract data
  otu_table <- get_otu_table(obj)
  tax_data <- get_tax_table(obj)

  # Detect spiked taxa
  if (!is.null(merged_spiked_species)) {
    # Check for matches in any taxonomic rank
    spiked_taxa <- unique(unlist(
      lapply(merged_spiked_species, function(taxon) {
        rownames(tax_data)[apply(tax_data, 1, function(row) any(row %in% taxon))]
      })
    ))
  } else if (!is.null(merged_spiked_hashcodes)) {
    spiked_taxa <- rownames(tax_data)[rownames(tax_data) %in% merged_spiked_hashcodes]
  }

  if (length(spiked_taxa) == 0) {
    stop("No matching spiked taxa found in the dataset.")
  }

  # Total reads
  total_reads <- data.frame(
    Sample = colnames(otu_table),
    Total_Reads = colSums(otu_table)
  )

  # Spiked reads
  spiked_reads <- data.frame(
    Sample = colnames(otu_table),
    Spiked_Reads = colSums(otu_table[spiked_taxa, , drop = FALSE])
  )

  # Merge and compute percentage
  merged_data <- dplyr::left_join(total_reads, spiked_reads, by = "Sample")
  merged_data <- dplyr::mutate(merged_data,
    Percentage = (Spiked_Reads / Total_Reads) * 100,
    Result = ifelse(Percentage >= passed_range[1] & Percentage <= passed_range[2],
      "passed", "failed"
    )
  )

  # Optional DOCX + CSV output
  if (!is.null(output_file)) {
    ft <- flextable::flextable(merged_data) |>
      flextable::fontsize(size = 10) |>
      flextable::font(part = "all", fontname = "Inconsolata") |>
      flextable::color(part = "header", color = "#6A0572") |>
      flextable::bold(part = "header") |>
      flextable::italic()

    flextable::save_as_docx(ft, path = output_file)
    utils::write.csv(merged_data, file = sub(".docx$", ".csv", output_file), row.names = FALSE)

    message("DOCX saved to: ", output_file)
    message("CSV saved to: ", sub(".docx$", ".csv", output_file))
  }

  return(merged_data)
}


# Example usage:
# Define the spiked species
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# calculate_spike_percentage(merged_TSE_sum, merged_spiked_species,
# passed_range = c(0.1, 20))


# merged_spiked_species<-"Dekkera_bruxellensis"
# result <- calculate_spike_percentage(merged_sum,
# merged_spiked_species,
# passed_range = c(0.1, 35))
# calculate_summary_stats_table(result)
# result$Percentage

# Define the spiked hashcodes
# merged_spiked_hashcodes <- c("hashcode1", "hashcode2")
# calculate_spike_percentage(physeq,
# merged_spiked_hashcodes = merged_spiked_hashcodes,
# passed_range = c(0.1, 10))
