#' @title Random Subsampling with Reduction Factor
#'
#' @description This function performs random subsampling on the OTU table of a phyloseq or
#' TreeSummarizedExperiment (TSE) object, reducing the counts of each ASV
#' (Amplicon Sequence Variant) by a specified reduction factor.
#' The resulting subsampled object is saved to a specified output file.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param reduction_factor A numeric value specifying the factor by which to reduce the counts. Default is 3.
#' @param output_file A character string specifying the output file name for the subsampled object. Default is "Less_subsampled_object.rds".
#' @return A subsampled object of the same class as the input.
#' @examples
#' \dontrun{
#' # Perform random subsampling on the phyloseq/TSE object with a reduction factor of 10
#' data("physeq_16SOTU", package="DspikeIn")
#' red <- random_subsample_WithReductionFactor(physeq_16SOTU, reduction_factor = 10)
#'
#' # Summarize the subsampled object
#' summary_stats <- summ_phyloseq_sampleID(red)
#' print(summary_stats)
#'
#' # TSE format
#' tse_16SOTU<- convert_phyloseq_to_tse("physeq_16SOTU")
#' red <- random_subsample_WithReductionFactor(tse_16SOTU, reduction_factor = 10)
#' summary_stats <- summ_phyloseq_sampleID(red)
#' }
#' @importFrom phyloseq otu_table phyloseq
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @export
random_subsample_WithReductionFactor <- function(obj, reduction_factor = 3, output_file = "Less_subsampled_object.rds") {
  suppressMessages({

    # Check object type and retrieve OTU table
    otu_table_mat <- get_otu_table(obj)

    if (is.null(otu_table_mat)) {
      stop("\U0000274C OTU table could not be extracted from the object.")
    }

    # Convert to data frame for manipulation
    otu_table_df <- as.data.frame(otu_table_mat)

    # Apply subsampling to each ASV in each sample
    for (i in seq_len(ncol(otu_table_df))) {
      asv_counts <- otu_table_df[, i]
      otu_table_df[, i] <- pmax(0, round(asv_counts / reduction_factor))
    }

    # Convert back to matrix
    otu_table_mod <- as.matrix(otu_table_df)

    # Reconstruct the object with the modified OTU table
    if (inherits(obj, "phyloseq")) {
      # Create a new phyloseq object
      subsampled_obj <- phyloseq::phyloseq(
        phyloseq::otu_table(otu_table_mod, taxa_are_rows = TRUE),
        phyloseq::sample_data(obj),
        phyloseq::tax_table(obj)
      )
    } else if (inherits(obj, "TreeSummarizedExperiment")) {
      # Create a new TSE object
      subsampled_obj <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = otu_table_mod),
        colData = SummarizedExperiment::colData(obj),
        rowData = SummarizedExperiment::rowData(obj),
        metadata = S4Vectors::metadata(obj)
      )
    } else {
      stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
    }

    # Save the subsampled object
    saveRDS(subsampled_obj, file = output_file)
    cat("Less_Subsampled object saved to:", output_file, "\n")
  })

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

