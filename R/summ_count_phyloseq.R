#' @title Summary Statistics of a Phyloseq or TSE Object
#' @description Computes overall summary statistics (mean, median, standard deviation, standard error, and quantiles)
#' for the OTU table in a `phyloseq` or `TreeSummarizedExperiment` (TSE) object.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing taxonomic and abundance data.
#' @return A data frame with overall summary statistics.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Summarize counts for the phyloseq object
#'   summary_stats_physeq <- summ_count_phyloseq(physeq_16SOTU)
#'
#'   # Convert phyloseq object to a TreeSummarizedExperiment (TSE)
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Summarize counts for the TSE object
#'   summary_stats_tse <- summ_count_phyloseq(tse_16SOTU)
#' }
#'
#' @importFrom matrixStats rowMeans2 rowMedians
#' @importFrom stats median sd quantile
#' @export
summ_count_phyloseq <- function(obj) {
  #  Convert TSE to phyloseq if needed
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  #  Extract OTU table using accessor
  otu_table <- get_otu_table(obj)

  #  Compute summary statistics
  overall_summary <- data.frame(
    Variable = c("Mean", "Median", "SD", "SE", "Q25", "Q50", "Q75"),
    Value = c(
      base::mean(matrixStats::rowMeans2(otu_table, na.rm = TRUE)),
      stats::median(matrixStats::rowMedians(otu_table, na.rm = TRUE)),
      stats::sd(matrixStats::rowMeans2(otu_table, na.rm = TRUE)),
      stats::sd(matrixStats::rowMeans2(otu_table, na.rm = TRUE)) / base::sqrt(nrow(otu_table)),
      stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.25),
      stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.50),
      stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.75)
    )
  )

  return(overall_summary)
}

# Example usage:
# summary_stats <- summ_count_phyloseq(physeq_ITSOTU)
# 
