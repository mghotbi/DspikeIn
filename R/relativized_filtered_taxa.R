#' @title Filter Taxa from a Phyloseq or TSE Object Based on Custom Thresholds
#'
#' @description This function filters taxa from a phyloseq or TreeSummarizedExperiment (TSE) object based on
#' custom thresholds for percentage of samples, mean abundance, count, and relative abundance.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param threshold_percentage A numeric value specifying the minimum percentage of samples in which
#' a taxon must be present to be retained. Default is 0.5.
#' @param threshold_mean_abundance A numeric value specifying the minimum mean abundance of a taxon
#' to be retained. Default is 0.001.
#' @param threshold_count A numeric value specifying the minimum count of a taxon in a sample
#' to be considered present. Default is 10.
#' @param threshold_relative_abundance A numeric value specifying the minimum relative abundance of a
#' taxon to be retained. Default is NULL.
#' @return A filtered phyloseq or TSE object containing only the taxa that meet the specified thresholds.
#' @importFrom phyloseq nsamples sample_sums filter_taxa
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Apply relative filtering on taxa
#'   FT <- relativized_filtered_taxa(
#'     physeq_16SOTU,
#'     threshold_percentage = 0.001,
#'     threshold_mean_abundance = 1,
#'     threshold_count = 5,
#'     threshold_relative_abundance = 0.001
#'   )
#' }
#' @export
relativized_filtered_taxa <- function(obj,
                                      threshold_percentage = 0.5,
                                      threshold_mean_abundance = 0.001,
                                      threshold_count = 10,
                                      threshold_relative_abundance = NULL) {
  # Get the OTU table
  otu_mat <- withCallingHandlers(
    get_otu_table(obj),
    message = function(m) invokeRestart("muffleMessage")
  )

  # Get the number of samples
  nsamples <- ncol(otu_mat)

  # Get sample sums
  sample_sum <- colSums(otu_mat)

  # Custom filter function
  filter_function <- function(x) {
    (sum(x > threshold_count) > nsamples * threshold_percentage) ||
      ((sum(x > threshold_count) > (nsamples * 0.1)) &&
        (mean(x / sample_sum) > threshold_mean_abundance) &&
        (!is.null(threshold_relative_abundance) &&
          max(x / sample_sum) > threshold_relative_abundance))
  }

  # Apply the filter
  taxa_to_keep <- apply(otu_mat, 1, filter_function)

  # Filter taxa and return the modified object
  if (inherits(obj, "phyloseq")) {
    return(prune_taxa(taxa_to_keep, obj))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(obj[taxa_to_keep, ])
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}


# Example usage with custom thresholds
# FT <- relativized_filtered_taxa(phyloseq_16SOTU, threshold_percentage = 0.001,
# threshold_mean_abundance = 1, threshold_count = 5, threshold_relative_abundance = 0.001)
