#' @title Tidy a Phyloseq or TreeSummarizedExperiment Object
#'
#' @description
#' Cleans and standardizes a microbiome object by:
#' - Removing taxonomic prefixes (e.g. `k__`, `p__`)
#' - Trimming whitespace
#' - Renaming `Domain` to `Kingdom` (if present)
#' - Removing taxa with "Chloroplast" or "Mitochondria" in **any** rank
#' - Filtering out taxa with zero total abundance
#' - Ensuring robust handling of both phyloseq and TreeSummarizedExperiment formats
#'
#' @param obj A `phyloseq::phyloseq` or `TreeSummarizedExperiment::TreeSummarizedExperiment` object.
#'
#' @return A cleaned and filtered object of the same class with updated taxonomy.
#'
#' @examples
#'
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   tidy_physeq <- tidy_phyloseq_tse(physeq_16SOTU)
#' }
#'
#' @importFrom phyloseq prune_taxa tax_table tax_table<-
#' @importFrom SummarizedExperiment rowData rowData<- assay
#' @importFrom methods is
#' @export
tidy_phyloseq_tse <- function(obj) {
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  tax_data <- get_tax_table(obj)

  if (is.null(tax_data) || ncol(tax_data) == 0) {
    warning("No taxonomy table found. Returning unmodified object.")
    return(obj)
  }

  # Rename Domain → Kingdom
  if ("Domain" %in% colnames(tax_data)) {
    colnames(tax_data)[colnames(tax_data) == "Domain"] <- "Kingdom"
  }

  # Clean taxonomy
  tax_data_clean <- as.data.frame(
    lapply(tax_data, function(col) {
      gsub("^[a-z]__\\s*", "", trimws(col))
    }),
    stringsAsFactors = FALSE
  )
  rownames(tax_data_clean) <- rownames(tax_data)

  # Remove unwanted taxa
  unwanted <- apply(tax_data_clean, 1, function(row) {
    any(grepl("chloroplast|mitochondria", row, ignore.case = TRUE))
  })

  # Get OTU table
  otu_matrix <- get_otu_table(obj)
  if (is.null(otu_matrix) || nrow(otu_matrix) == 0) {
    warning("No OTU table found. Returning unmodified object.")
    return(obj)
  }

  non_zero <- rowSums(otu_matrix) > 0
  keep_taxa_ids <- rownames(otu_matrix)[non_zero & !unwanted]

  # --- Subset and synchronize
  if (inherits(obj, "phyloseq")) {
    obj <- phyloseq::prune_taxa(keep_taxa_ids, obj)

    # Match cleaned taxonomy to remaining taxa
    tax_clean_subset <- tax_data_clean[keep_taxa_ids, , drop = FALSE]
    stopifnot(identical(keep_taxa_ids, rownames(tax_clean_subset)))

    phyloseq::tax_table(obj) <- phyloseq::tax_table(as.matrix(tax_clean_subset))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- obj[keep_taxa_ids, ]
    SummarizedExperiment::rowData(obj) <- tax_data_clean[keep_taxa_ids, , drop = FALSE]
  }

  return(obj)
}


# Usage Example
# tidy_physeq <- tidy_phyloseq_tse(M18)
# tidy_Mt_tse <- tidy_phyloseq_tse(mt)
# taxonomy_table <- SummarizedExperiment::rowData(tidy_M19_tse)
# phylo_tree <- S4Vectors::metadata(tidy_M19_tse)$tree
