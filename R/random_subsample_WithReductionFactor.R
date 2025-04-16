#' @title Random Subsampling with Reduction Factor
#'
#' @description
#' Performs random subsampling on the OTU table of a `phyloseq` or
#' `TreeSummarizedExperiment` (TSE) object by dividing each ASV count
#' by a specified reduction factor and rounding down to the nearest whole number.
#' Optionally, the result can be saved to disk as an `.rds` file.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param reduction_factor A numeric value \eqn{\ge 1} to reduce counts. Default is 3.
#' @param output_file Optional. A character string specifying the `.rds` file path
#'   to save the result. If `NULL`, no file will be saved. Default is `NULL`.
#'
#' @return A subsampled object of the same class as the input (`phyloseq` or `TreeSummarizedExperiment`).
#'
#' @examples
#' data("physeq_16SOTU", package = "DspikeIn")
#' red <- random_subsample_WithReductionFactor(physeq_16SOTU, reduction_factor = 10)
#' summary_stats <- summ_phyloseq_sampleID(red)
#' print(summary_stats)
#'
#' # Subsampling TSE object
#' tse_16SOTU <- convert_phyloseq_to_tse("physeq_16SOTU")
#' red <- random_subsample_WithReductionFactor(tse_16SOTU, reduction_factor = 10)
#' summary_stats <- summ_phyloseq_sampleID(red)
#' print(summary_stats)
#'
#' @importFrom phyloseq otu_table phyloseq sample_data tax_table
#' @importFrom SummarizedExperiment assay SummarizedExperiment colData rowData
#' @importFrom S4Vectors metadata
#' @export
random_subsample_WithReductionFactor <- function(obj,
                                                 reduction_factor = 3,
                                                 output_file = NULL) {
  # Validate reduction_factor
  if (!is.numeric(reduction_factor) || reduction_factor <= 0) {
    stop("`reduction_factor` must be a positive numeric value.")
  }

  # Extract OTU/count matrix
  otu_table_mat <- get_otu_table(obj)

  if (is.null(otu_table_mat)) {
    stop("OTU table could not be extracted from the input object.")
  }

  # Convert to data frame for manipulation
  otu_table_df <- as.data.frame(otu_table_mat)

  # Subsample by dividing and rounding counts
  for (i in seq_len(ncol(otu_table_df))) {
    asv_counts <- otu_table_df[, i]
    otu_table_df[, i] <- pmax(0, round(asv_counts / reduction_factor))
  }

  # Convert back to matrix
  otu_table_mod <- as.matrix(otu_table_df)

  # Rebuild the object depending on class
  if (inherits(obj, "phyloseq")) {
    subsampled_obj <- phyloseq::phyloseq(
      phyloseq::otu_table(otu_table_mod, taxa_are_rows = TRUE),
      phyloseq::sample_data(obj),
      phyloseq::tax_table(obj)
    )
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    subsampled_obj <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = otu_table_mod),
      colData = SummarizedExperiment::colData(obj),
      rowData = SummarizedExperiment::rowData(obj),
      metadata = S4Vectors::metadata(obj)
    )
  } else {
    stop("Unsupported object type. Must be `phyloseq` or `TreeSummarizedExperiment`.")
  }

  # Optionally save output
  if (!is.null(output_file)) {
    if (!is.character(output_file) || length(output_file) != 1) {
      stop("`output_file` must be a single string or NULL.")
    }

    saveRDS(subsampled_obj, file = output_file)
    message("Subsampled object saved to: ", output_file)
  }

  return(subsampled_obj)
}

# Example usage:
# Perform random subsampling with a reduction factor of 10
# red <- random_subsample_WithReductionFactor(spiked_16S, reduction_factor = 10)
# Summarize the subsampled phyloseq object
# summ_phyloseq_sampleID(red)
# red_physeq <- random_subsample_WithReductionFactor(TSE_obj, reduction_factor = 10)
# summary_stats <- summ_phyloseq_sampleID(red_physeq)
# print(summary_stats)
