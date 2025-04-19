#' @title Adjust Prevalence in a Microbiome Object
#'
#' @description
#' Removes low-prevalence taxa from a `phyloseq` or `TreeSummarizedExperiment` object based
#' on user-specified prevalence thresholds derived from taxa abundance statistics.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param method Character. Threshold method: one of \code{"min"}, \code{"mean"}, \code{"median"}, or \code{"max"}.
#' @param output_file Optional. Character. Path to save the adjusted object as `.rds`.
#' Default is `NULL`, meaning no file will be saved unless explicitly provided.
#'
#' @return The adjusted object of the same class as `obj`.
#'
#' @importFrom phyloseq prune_taxa taxa_sums sample_sums otu_table
#' @importFrom SummarizedExperiment assay
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   ## Adjust prevalence in phyloseq
#'   adjusted_physeq <- adjusted_prevalence(physeq_16SOTU, method = "mean")
#'
#'   ## Convert to TreeSummarizedExperiment
#'   tse_obj <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   ## Adjust prevalence in TSE
#'   adjusted_tse <- adjusted_prevalence(tse_obj, method = "median")
#' }
#'
#' @export
adjusted_prevalence <- function(obj, method = "min", output_file = NULL) {
  # Validate method
  method <- match.arg(method, c("min", "mean", "median", "max"))

  # Extract data depending on object type
  if (inherits(obj, "phyloseq")) {
    otu_mat <- as.matrix(phyloseq::otu_table(obj))
    taxa_sums_vec <- phyloseq::taxa_sums(obj)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    otu_mat <- as.matrix(SummarizedExperiment::assay(obj))
    taxa_sums_vec <- rowSums(otu_mat)
  } else {
    stop("Object must be of class 'phyloseq' or 'TreeSummarizedExperiment'.")
  }

  # Calculate abundance threshold
  abundance_threshold <- switch(method,
    "min" = min(taxa_sums_vec),
    "mean" = mean(taxa_sums_vec),
    "median" = median(taxa_sums_vec),
    "max" = max(taxa_sums_vec)
  )

  # Calculate prevalence threshold (default = present in approx. 10% of samples)

  prevalence_counts <- rowSums(otu_mat > 0)
  prevalence_threshold <- ceiling(0.1 * ncol(otu_mat))
  # (default = present in >= 10% of samples)
  message("Prevalence threshold = taxa detected in >= ", prevalence_threshold, " samples.")
  message("Abundance threshold (", method, ") = ", round(abundance_threshold, 2), " reads.")

  # Filter
  selected <- which((prevalence_counts >= prevalence_threshold) & (taxa_sums_vec >= abundance_threshold))

  if (length(selected) == 0) stop("No taxa passed the chosen thresholds.")

  message("Number of retained taxa: ", length(selected))

  # Apply filtering
  if (inherits(obj, "phyloseq")) {
    obj_adj <- phyloseq::prune_taxa(rownames(otu_mat)[selected], obj)
  } else {
    obj_adj <- obj[selected, ]
  }

  # Optional output
  if (!is.null(output_file)) {
    saveRDS(obj_adj, file = output_file)
    message("Adjusted object saved to: ", output_file)
  }

  return(obj_adj)
}


# Usage Example:
# adjusted_physeq <- adjusted_prevalence(physeq_16SOTU, method = "min")
# adjusted_physeq <- adjusted_prevalence(physeq_ITSOTU, method = "min")
# tse_16SOTU<-convert_phyloseq_to_tse(physeq_16SOTU)
# adjusted_tse <- adjusted_prevalence(tse_16SOTU, method = "min")
