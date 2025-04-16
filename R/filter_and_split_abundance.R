#' @title Filter and Split Abundance Data by Threshold
#' @description
#' Filters low-abundance taxa from a `phyloseq` or `TreeSummarizedExperiment` object,
#' splits data into high- and low-abundance groups, and optionally saves results.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param threshold A numeric value indicating the mean abundance threshold for filtering.
#' @param output_prefix A character string specifying the filename prefix for saving.
#'                      If `NULL`, files will not be saved. Default is `NULL`.
#'
#' @return A named list with components:
#' \describe{
#'   \item{high}{Subset with taxa having mean abundance > threshold}
#'   \item{low}{Subset with taxa having mean abundance <= threshold}
#' }
#' Each returned element contains a single object of the same class as the input (either \code{phyloseq} or \code{TreeSummarizedExperiment}).
#'
#' @examples
#' data("physeq_ITSOTU", package = "DspikeIn")
#'
#' # Return results without saving
#' output <- filter_and_split_abundance(physeq_ITSOTU, threshold = 0.05)
#'
#' # With saving
#' tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
#'
#' output <- filter_and_split_abundance(tse_ITSOTU,
#'   threshold = 0.05,
#'   output_prefix = file.path(tempdir(), "abund")
#' )
#'
#' @importFrom phyloseq prune_taxa
#' @export
filter_and_split_abundance <- function(obj, threshold = 0.01, output_prefix = NULL) {
  # Validate input
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  }

  is_physeq <- inherits(obj, "phyloseq")

  # Get abundance matrix
  otu <- get_otu_table(obj)

  # Identify high and low abundance taxa
  high_otu <- otu[rowMeans(otu) > threshold, , drop = FALSE]
  low_otu <- otu[rowMeans(otu) <= threshold, , drop = FALSE]

  if (is_physeq) {
    high_obj <- phyloseq::prune_taxa(rownames(high_otu), obj)
    low_obj <- phyloseq::prune_taxa(rownames(low_otu), obj)
  } else {
    high_obj <- obj[rownames(high_otu), ]
    low_obj <- obj[rownames(low_otu), ]
  }

  # Save to files if prefix provided
  if (!is.null(output_prefix)) {
    saveRDS(high_obj, paste0(output_prefix, "_high.rds"))
    saveRDS(low_obj, paste0(output_prefix, "_low.rds"))

    message("\n[OK] High- and low-abundance objects saved:")
    message("High: ", output_prefix, "_high.rds")
    message("Low : ", output_prefix, "_low.rds")
  }

  return(list(
    high = high_obj,
    low = low_obj
  ))
}
