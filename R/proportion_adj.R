#' @title Proportionally Adjust Abundance
#' @description This function normalizes the abundance data in a `phyloseq` or `TreeSummarizedExperiment`
#' object by adjusting each sample's counts based on a total value, typically the maximum total
#' sequence count across all samples. The adjusted counts are then rounded to the nearest integer.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param output_file A character string specifying the output file name for the adjusted object.
#' If `NULL`, the object is not saved. Default is "proportion_adjusted.rds".
#' @return A modified object of the same class (`phyloseq` or `TreeSummarizedExperiment`) with proportionally adjusted and rounded abundance data.
#'
#' @details
#' This function extracts the OTU table (or assay in `TSE`), normalizes it based on the sample sums,
#' and updates the original object while maintaining its structure.
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load phyloseq object
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   normalized_physeq <- proportion_adj(
#'     physeq_16SOTU,
#'     output_file = file.path(tempdir(), "proportion_adjusted_physeq.rds")
#'   )
#'   print(normalized_physeq)
#'
#'   # Convert to TSE and apply
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   normalized_tse <- proportion_adj(
#'     tse_16SOTU,
#'     output_file = file.path(tempdir(), "proportion_adjusted_tse.rds")
#'   )
#'   print(normalized_tse)
#' }
#' @importFrom phyloseq otu_table<-
#' @importFrom SummarizedExperiment assay<-
#' @importFrom phyloseq sample_sums transform_sample_counts
#' @importFrom S4Vectors metadata
#' @export
proportion_adj <- function(obj, output_file = "proportion_adjusted.rds") {
  message("Starting proportional adjustment...")

  # Extract OTU table using accessor function
  otu_matrix <- get_otu_table(obj)
  if (is.null(otu_matrix)) stop("Error: OTU table is missing.")

  # Extract sample sums using accessor
  sample_totals <- get_sample_sums(obj)
  if (is.null(sample_totals)) stop("Error: Sample sums could not be calculated.")

  # Compute normalization function
  normf <- function(x, tot = max(sample_totals)) {
    tot * x / sum(x)
  }

  message("Applying normalization based on maximum total sequence count...")

  # Apply normalization to sample counts/abundance
  if (inherits(obj, "phyloseq")) {
    obj <- phyloseq::transform_sample_counts(obj, normf)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    otu_matrix <- apply(otu_matrix, 2, normf)
  }

  # Round the abundance counts
  message("Rounding OTU table values...")
  otu_matrix <- round(otu_matrix, digits = 0)

  # Update OTU table in the object
  if (inherits(obj, "phyloseq")) {
    phyloseq::otu_table(obj) <- phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::assay(obj) <- otu_matrix
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  # Save if output file is provided
  if (!is.null(output_file)) {
    message("Saving adjusted object to: ", output_file)
    saveRDS(obj, file = output_file)
  }

  return(obj)
}

#' @title Extract OTU Table from Object
#' @description Retrieves the OTU table from a `phyloseq` or `TreeSummarizedExperiment` object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A matrix containing OTU count data.
#' @importFrom phyloseq otu_table
#' @importFrom SummarizedExperiment assay
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

#' @title Extract Sample Sums from Object
#' @description Retrieves the total sequence counts per sample from a `phyloseq` or `TreeSummarizedExperiment` object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A numeric vector of total counts per sample.
#' @importFrom phyloseq sample_sums
#' @export
get_sample_sums <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    return(phyloseq::sample_sums(obj))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(colSums(SummarizedExperiment::assay(obj)))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}
