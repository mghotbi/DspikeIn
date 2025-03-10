#' @title Tidy a Phyloseq or TreeSummarizedExperiment Object
#' @description Cleans and standardizes a microbiome dataset, supporting both
#' `phyloseq` and `TreeSummarizedExperiment`. Performs:
#' - Standardization of taxonomic ranks (if available)
#' - Removal of leading/trailing whitespace in taxa names
#' - Filtering out zero-count taxa
#' - Exclusion of "Chloroplast" and "Mitochondria" classifications (if applicable)
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A cleaned and tidied object of the same class.
#'
#' @details
#' This function standardizes taxonomic ranks, removes unnecessary whitespace,
#' and filters unwanted classifications, ensuring consistency for downstream analysis.
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Convert phyloseq object to tidy format
#'   tidy_physeq <- tidy_phyloseq_tse(physeq_16SOTU)
#'
#'   # TreeSummarizedExperiment (TSE) object
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Convert TSE object to tidy format
#'   tidy_tse <- tidy_phyloseq_tse(tse_16SOTU)
#' }
#' }
#' @importFrom phyloseq prune_taxa
#' @importFrom SummarizedExperiment assay
#' @export
tidy_phyloseq_tse <- function(obj) {
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  # Extract taxonomy table
  tax_data <- get_tax_table(obj)

  if (is.null(tax_data) || ncol(tax_data) == 0) {
    warning("\U000026A0 No taxonomy table found. Returning unmodified object.")
    return(obj)
  }

  # Fix taxa names by removing "__" prefixes (e.g., "k__Bacteria" → "Bacteria")
  tax_data <- as.data.frame(lapply(tax_data, function(col) gsub("[a-z]__\\s*", "", col)), stringsAsFactors = FALSE)

  # Define standard taxonomic ranks
  required_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  available_ranks <- intersect(required_ranks, colnames(tax_data))

  # Ensure taxonomic table maintains structure
  if (length(available_ranks) > 0) {
    tax_data <- tax_data[, available_ranks, drop = FALSE]
  } else {
    warning("\U000026A0 No standard taxonomic ranks found. Proceeding with available data.")
  }

  # Trim whitespace from taxonomic names
  tax_data <- as.data.frame(lapply(tax_data, trimws), stringsAsFactors = FALSE)

  # Extract OTU table
  otu_matrix <- get_otu_table(obj)

  if (is.null(otu_matrix) || nrow(otu_matrix) == 0) {
    warning("\U000026A0 No OTU table found. Returning unmodified object.")
    return(obj)
  }

  # Identify and remove zero-count taxa
  non_zero_taxa <- rowSums(otu_matrix) > 0

  # Remove unwanted taxa based on taxonomy
  if ("Class" %in% colnames(tax_data)) {
    non_zero_taxa <- non_zero_taxa & tax_data$Class != "Chloroplast"
  }
  if ("Family" %in% colnames(tax_data)) {
    non_zero_taxa <- non_zero_taxa & tax_data$Family != "Mitochondria"
  }

  # Apply filtering to the object
  if (inherits(obj, "phyloseq")) {
    obj <- phyloseq::prune_taxa(rownames(otu_matrix)[non_zero_taxa], obj)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- obj[rownames(otu_matrix)[non_zero_taxa], ]
  }

  return(obj)
}

#Usage Example
# tidy_physeq <- tidy_phyloseq_tse(M18)
# tidy_Mt_tse <- tidy_phyloseq_tse(mt)
#taxonomy_table <- SummarizedExperiment::rowData(tidy_M19_tse)
# phylo_tree <- S4Vectors::metadata(tidy_M19_tse)$tree


