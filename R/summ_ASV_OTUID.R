#' @title Summarize ASV Data Based on ASV_ID
#' @description This function generates summary statistics (mean, median, standard deviation,
#' standard error, and quantiles) for each ASV (Amplicon Sequence Variant)
#' in a `phyloseq` or `TreeSummarizedExperiment` object.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing ASV abundance data.
#' @return A data frame containing summary statistics for each ASV.
#'
#' @details
#' This function extracts the OTU table (or assay in `TSE`), computes per-ASV statistics,
#' and returns a tidy summary data frame.
#'
#' @examples
#' # Example with a phyloseq object
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_ITSOTU", package = "DspikeIn")
#'   summary_physeq <- summ_ASV_OTUID(physeq_ITSOTU)
#'
#'   # Example with a TreeSummarizedExperiment object
#'   tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
#'   summary_tse <- summ_ASV_OTUID(tse_ITSOTU)
#' }
#'
#' @importFrom stats quantile sd
#' @importFrom phyloseq otu_table
#' @importFrom SummarizedExperiment assay
#' @export
summ_ASV_OTUID <- function(obj) {
  # Extract OTU table using accessor function
  otu_table <- get_otu_table(obj)
  if (is.null(otu_table)) stop("Error: OTU table is missing.")

  # Calculate summary statistics for each ASV
  summary_stats <- data.frame(
    ASV_ID = rownames(otu_table), # ASV identifier
    Mean = rowMeans(otu_table, na.rm = TRUE), # Mean abundance
    Median = apply(otu_table, 1, median, na.rm = TRUE), # Median abundance
    SD = apply(otu_table, 1, stats::sd, na.rm = TRUE), # Standard deviation
    SE = apply(otu_table, 1, function(x) stats::sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))), # Standard error
    Q25 = apply(otu_table, 1, stats::quantile, probs = 0.25, na.rm = TRUE), # 25th percentile (Q1)
    Q50 = apply(otu_table, 1, stats::quantile, probs = 0.5, na.rm = TRUE), # 50th percentile (Q2 / median)
    Q75 = apply(otu_table, 1, stats::quantile, probs = 0.75, na.rm = TRUE) # 75th percentile (Q3)
  )

  return(summary_stats)
}

#' @title Extract OTU Table from Object
#' @description Retrieves the OTU table from a `phyloseq` or `TreeSummarizedExperiment` object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A matrix containing OTU count data.
#' @export
get_otu_table <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    return(as.matrix(phyloseq::otu_table(obj)))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(as.matrix(SummarizedExperiment::assay(obj)))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}
# 
