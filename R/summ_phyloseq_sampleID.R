#' @title Generate Summary Statistics for Each Sample
#'
#' @description Calculates summary statistics (mean, median, standard deviation, standard error, and quartiles)
#' for each sample in a \code{phyloseq} or \code{TreeSummarizedExperiment} object.
#'
#' @param obj A \code{phyloseq} or \code{TreeSummarizedExperiment} object containing microbiome data.
#'
#' @return A data frame containing summary statistics per sample, with columns:
#' \itemize{
#'   \item Sample_ID
#'   \item Mean
#'   \item Median
#'   \item Standard Deviation (SD)
#'   \item Standard Error (SE)
#'   \item Q25 (25th percentile)
#'   \item Q50 (Median)
#'   \item Q75 (75th percentile)
#' }
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Summarize the phyloseq object
#'   summary_stats_physeq <- summ_phyloseq_sampleID(physeq_16SOTU)
#'   print(summary_stats_physeq)
#'
#'   # Convert to TreeSummarizedExperiment
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   summary_stats_tse <- summ_phyloseq_sampleID(tse_16SOTU)
#'   print(summary_stats_tse)
#' }
#'
#' @importFrom stats quantile sd
#' @export
summ_phyloseq_sampleID <- function(obj) {
  # Extract OTU table (as matrix)
  otu_matrix <- get_otu_table(obj)

  # Validate input
  if (!is.matrix(otu_matrix) || ncol(otu_matrix) == 0) {
    stop("Error: OTU table is missing or has no samples.")
  }

  # Calculate summary statistics
  summary_stats <- data.frame(
    Sample_ID = colnames(otu_matrix),
    Mean = apply(otu_matrix, 2, mean, na.rm = TRUE),
    Median = apply(otu_matrix, 2, median, na.rm = TRUE),
    SD = apply(otu_matrix, 2, sd, na.rm = TRUE),
    SE = apply(otu_matrix, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    Q25 = apply(otu_matrix, 2, quantile, probs = 0.25, na.rm = TRUE),
    Q50 = apply(otu_matrix, 2, quantile, probs = 0.50, na.rm = TRUE),
    Q75 = apply(otu_matrix, 2, quantile, probs = 0.75, na.rm = TRUE)
  )

  return(summary_stats)
}

# Example usage:
# Generate summary statistics for a phyloseq object
# summary_stats <- summ_phyloseq_sampleID(phyloseq-obj)
# summary_stats <- summ_phyloseq_sampleID(tse-obj)
