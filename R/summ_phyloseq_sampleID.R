#' @title Generate Summary Statistics for Each Sample
#' @description Calculates summary statistics (mean, median, standard deviation, standard error, and quartiles)
#'              for each sample in a `phyloseq` or `TreeSummarizedExperiment` object.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @return A data frame containing summary statistics per sample, with columns:
#'         Sample_ID, Mean, Median, Standard Deviation, Standard Error, Q25, Q50, Q75.
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Summarize the phyloseq object by sample ID
#'   summary_stats_physeq <- summ_phyloseq_sampleID(physeq_16SOTU)
#'   print(summary_stats_physeq)
#'
#'   # Convert phyloseq object to a TreeSummarizedExperiment (TSE)
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Summarize the TSE object by sample ID
#'   summary_stats_tse <- summ_phyloseq_sampleID(tse_16SOTU)
#'   print(summary_stats_tse)
#' }
#' }
#'
#' @importFrom stats quantile sd
#' @export
summ_phyloseq_sampleID <- function(obj) {
  suppressMessages({

    #  Extract OTU table using accessor function
    otu_matrix <- get_otu_table(obj)

    #  Ensure data is in the correct orientation (samples in columns)
    if (ncol(otu_matrix) == 0) {
      stop("OTU table is empty or incorrectly formatted.")
    }

    # Compute summary statistics for each sample
    summary_stats <- data.frame(
      Sample_ID = colnames(otu_matrix),
      Mean = apply(otu_matrix, 2, mean, na.rm = TRUE),
      Median = apply(otu_matrix, 2, median, na.rm = TRUE),
      SD = apply(otu_matrix, 2, stats::sd, na.rm = TRUE),
      SE = apply(otu_matrix, 2, function(x) stats::sd(x, na.rm = TRUE) / sqrt(length(x))),
      Q25 = apply(otu_matrix, 2, stats::quantile, probs = 0.25, na.rm = TRUE),
      Q50 = apply(otu_matrix, 2, stats::quantile, probs = 0.50, na.rm = TRUE),
      Q75 = apply(otu_matrix, 2, stats::quantile, probs = 0.75, na.rm = TRUE)
    )

    return(summary_stats)
  })
}

# Example usage:
# Generate summary statistics for a phyloseq object
# summary_stats <- summ_phyloseq_sampleID(phyloseq-obj)
# summary_stats <- summ_phyloseq_sampleID(tse-obj)
