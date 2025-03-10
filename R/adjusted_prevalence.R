#' @title Adjust Prevalence in a Microbiome Object
#' @description Adjusts the prevalence of ASVs based on a specified method
#'   and prunes those that do not meet the prevalence threshold.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param method A character string specifying the method ("min", "mean", "median", "max").
#' @param output_file A character string specifying the file path to save the adjusted object.
#' @return An adjusted object with pruned taxa.
#' @importFrom phyloseq taxa_sums prune_taxa sample_sums
#' @importFrom SummarizedExperiment assay colData
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Adjust prevalence using the "min" method
#'   adjusted_physeq <- adjusted_prevalence(physeq_16SOTU, method = "min")
#'   print(adjusted_physeq)
#'
#'   # Convert phyloseq to TreeSummarizedExperiment (TSE)
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Adjust prevalence using the "max" method on a TSE object
#'   adjusted_tse <- adjusted_prevalence(tse_16SOTU, method = "max")
#'   print(adjusted_tse)
#' }
#' }
#' @export
adjusted_prevalence <- function(obj, method = "min", output_file = "adjusted_prevalence.rds") {
  suppressMessages({
    # Ensure method is valid
    method <- tolower(method)
    valid_methods <- c("min", "mean", "median", "max")
    if (!method %in% valid_methods) {
      stop("\U0000274C Invalid method. Choose from 'min', 'mean', 'median', or 'max'.")
    }

    # Extract abundance data
    if (inherits(obj, "phyloseq")) {
      reads <- switch(method,
                      "min" = min(phyloseq::taxa_sums(obj)),
                      "mean" = mean(phyloseq::taxa_sums(obj)),
                      "median" = median(phyloseq::taxa_sums(obj)),
                      "max" = max(phyloseq::taxa_sums(obj)))
      sample_sums <- phyloseq::sample_sums(obj)
      otu_matrix <- as.matrix(phyloseq::otu_table(obj))
    } else if (inherits(obj, "TreeSummarizedExperiment")) {
      abundance_data <- as.matrix(SummarizedExperiment::assay(obj))
      reads <- switch(method,
                      "min" = min(rowSums(abundance_data)),
                      "mean" = mean(rowSums(abundance_data)),
                      "median" = median(rowSums(abundance_data)),
                      "max" = max(rowSums(abundance_data)))
      sample_sums <- colSums(abundance_data)
      otu_matrix <- abundance_data
    } else {
      stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
    }

    # Print chosen method and reads
    cat("Number of reads chosen by method", method, ":", reads, "\n")

    # Compute prevalence
    prevalence_counts <- rowSums(otu_matrix > 0)
    otu_df <- data.frame(prev = prevalence_counts, sums = rowSums(otu_matrix))
    otu_df <- otu_df[order(-otu_df$prev, -otu_df$sums), ]

    # Define threshold and number of OTUs to keep
    prevalence_threshold <- 0.1 * reads
    nOTUs <- sum(otu_df$prev >= prevalence_threshold)
    cat("Prevalence threshold:", prevalence_threshold, "\n")
    cat("Number of OTUs to keep:", nOTUs, "\n")

    # Prune low-prevalence taxa
    selected_otus <- rownames(otu_df)[1:nOTUs]

    if (inherits(obj, "phyloseq")) {
      obj_adj <- phyloseq::prune_taxa(selected_otus, obj)
    } else if (inherits(obj, "TreeSummarizedExperiment")) {
      obj_adj <- obj[row.names(otu_matrix) %in% selected_otus, ]
    }

    # Save and return adjusted object
    saveRDS(obj_adj, file = output_file)
    cat("Adjusted object saved to:", output_file, "\n")

    return(obj_adj)
  })
}

#Usage Example:
# adjusted_physeq <- adjusted_prevalence(physeq_16SOTU, method = "min")
# adjusted_physeq <- adjusted_prevalence(physeq_ITSOTU, method = "min")
#tse_16SOTU<-convert_phyloseq_to_tse(physeq_16SOTU)
#adjusted_physeq <- adjusted_prevalence(tse_16SOTU, method = "min")

