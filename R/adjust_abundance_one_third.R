#' @title Adjust Abundance by One-Third
#' @description This function normalizes the abundance data in a `phyloseq` or `TreeSummarizedExperiment`
#' object by dividing each value by a specified factor.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param factor A numeric value specifying the factor by which to divide the abundance data. Default is 3.
#' @param output_file A character string specifying the output file name for the adjusted object. Default is NULL.
#' @return A modified object of the same class (`phyloseq` or `TreeSummarizedExperiment`) with adjusted abundances.
#'
#' @details
#' This function extracts the OTU table (or assay in `TSE`), normalizes it by dividing by `factor`,
#' and updates the original object while maintaining its structure.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   adjusted_physeq <- adjust_abundance_one_third(physeq_16SOTU, factor = 3)
#'
#'   # Example with a TreeSummarizedExperiment object
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   adjusted_tse <- adjust_abundance_one_third(tse_16SOTU, factor = 3)
#' }
#' }
#'
#' @importFrom phyloseq otu_table<-
#' @importFrom SummarizedExperiment assay<-
#' @importFrom S4Vectors metadata
#' @export
adjust_abundance_one_third <- function(obj, factor = 3, output_file = NULL) {
  message("Starting abundance adjustment...")

  # Extract OTU table
  otu_matrix <- get_otu_table(obj)
  if (is.null(otu_matrix)) stop("Error: OTU table is missing.")

  message("Dividing OTU table by factor ", factor, "...")
  otu_matrix <- otu_matrix / factor  # Normalize abundance

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
    message("\U0001F4BE Saving adjusted object to: ", output_file)
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
    stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}
