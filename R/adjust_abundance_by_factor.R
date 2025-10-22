#' @title Adjust Abundance by a Custom Factor
#'
#' @description
#' This function normalizes microbial abundance data in a `phyloseq` or `TreeSummarizedExperiment`
#' object by dividing all abundance values by a user-defined factor.
#' It supports both data structures and preserves their internal metadata.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome count data.
#' @param factor A numeric value \eqn{> 0} specifying the divisor applied to abundance counts.
#'   Default is 3 (historically used as the one-third normalization factor).
#' @param output_file A character string specifying a path to save the adjusted object as `.rds`.
#'   Default is `NULL` (no file written).
#'
#' @return
#' An adjusted object of the same class (`phyloseq` or `TreeSummarizedExperiment`),
#' where the abundance values are divided by the specified `factor`.
#'
#' @details
#' This function extracts the OTU table (for `phyloseq`) or assay matrix (for `TSE`),
#' divides all abundance values by the provided `factor`, and reinserts the adjusted
#' table while maintaining the full metadata structure.
#'
#' Historically, the function name `adjust_abundance_one_third()` referred to a fixed
#' normalization by 3. This generalized version preserves the same function name for
#' backward compatibility and allows any user-defined divisor.
#'
#' For clarity, an alias function `adjust_abundance_by_factor()` is provided.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Adjust phyloseq object
#'   adjusted_physeq <- adjust_abundance_one_third(physeq_16SOTU, factor = 3)
#'
#'   # Convert to TSE and adjust
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   adjusted_tse <- adjust_abundance_one_third(tse_16SOTU, factor = 2)
#'
#'   # Using the alias
#'   adjusted_physeq2 <- adjust_abundance_by_factor(physeq_16SOTU, factor = 5)
#' }
#'
#' @seealso
#' \code{\link{convert_phyloseq_to_tse}}
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @export
adjust_abundance_one_third <- function(obj, factor = 3, output_file = NULL) {
  message("Starting abundance adjustment...")
  
  # -------------------------------
  # Validate inputs
  # -------------------------------
  if (!is.numeric(factor) || factor <= 0)
    stop("'factor' must be numeric and greater than zero.")
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment")))
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  
  # -------------------------------
  # Extract abundance table
  # -------------------------------
  otu_matrix <- get_otu_table(obj)
  if (is.null(otu_matrix))
    stop("No abundance table found in the object.")
  
  # -------------------------------
  # Perform adjustment
  # -------------------------------
  message("Dividing all abundance values by factor = ", factor, " ...")
  otu_matrix <- otu_matrix / factor
  
  # -------------------------------
  # Update object
  # -------------------------------
  if (inherits(obj, "phyloseq")) {
    are_rows <- phyloseq::taxa_are_rows(obj)
    phyloseq::otu_table(obj) <- phyloseq::otu_table(
      otu_matrix,
      taxa_are_rows = are_rows
    )
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::assay(obj) <- otu_matrix
  }
  
  # -------------------------------
  # Save if requested
  # -------------------------------
  if (!is.null(output_file)) {
    message("Saving adjusted object to: ", output_file)
    saveRDS(obj, file = output_file)
  }
  
  message("Abundance adjustment complete.")
  invisible(obj)
}


#' @rdname adjust_abundance_one_third
#' @export
adjust_abundance_by_factor <- function(obj, factor = 3, output_file = NULL) {
  adjust_abundance_one_third(obj, factor = factor, output_file = output_file)
}


#' @title Extract OTU Table from Object
#' @description
#' Retrieves the OTU table or assay matrix from a `phyloseq` or `TreeSummarizedExperiment` object.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#'
#' @return A numeric matrix containing OTU count data.
#'
#' @details
#' Ensures consistent extraction across object classes and automatically transposes
#' the `phyloseq` matrix if taxa are not stored as rows.
#'
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom SummarizedExperiment assay
#' @export
get_otu_table <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    mat <- as.matrix(phyloseq::otu_table(obj))
    if (!phyloseq::taxa_are_rows(obj)) mat <- t(mat)
    return(mat)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(as.matrix(SummarizedExperiment::assay(obj)))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}
