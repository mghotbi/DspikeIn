#' @title Convert a `phyloseq` Object to a `TreeSummarizedExperiment`
#' @description Converts a `phyloseq` object into a `TreeSummarizedExperiment`
#'              (TSE) object, ensuring phylogenetic tree and reference sequences
#'              are preserved in `rowTree` and `metadata` if they exist.
#'
#' @param physeq A `phyloseq` object containing:
#'   \itemize{
#'     \item **OTU Table**: Feature abundance matrix.
#'     \item **Taxonomy Table** (optional): Taxonomic classification of features.
#'     \item **Sample Metadata** (optional): Experimental sample data.
#'     \item **Phylogenetic Tree** (optional).
#'     \item **Reference Sequences** (optional).
#'   }
#' @return A `TreeSummarizedExperiment` object preserving all components.
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'    data("physeq_16SOTU", package = "DspikeIn")
#'
#'    # Convert phyloseq object to TreeSummarizedExperiment
#'    tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#' }
#' }
#' @importFrom phyloseq otu_table tax_table sample_data phy_tree refseq taxa_are_rows sample_names
#' @importFrom microbiome meta
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom ape as.phylo is.rooted
#' @importFrom S4Vectors metadata
#' @export
convert_phyloseq_to_tse <- function(physeq) {
  if (!inherits(physeq, "phyloseq")) {
    stop(" Error: Input must be a valid 'phyloseq' object.")
  }

  # --- Extract OTU Table (counts) ---
  otu_matrix <- as(phyloseq::otu_table(physeq), "matrix")

  # Ensure taxa are in rows
  if (!phyloseq::taxa_are_rows(physeq)) {
    otu_matrix <- t(otu_matrix)
  }

  # --- Extract Taxonomy Table (if available) ---
  tax_data <- tryCatch(
    as.data.frame(phyloseq::tax_table(physeq)),
    error = function(e) NULL
  )

  # If taxonomy is missing, create an empty data frame with at least one column
  if (is.null(tax_data)) {
    tax_data <- data.frame(Dummy_Taxonomy = rep(NA, nrow(otu_matrix)), row.names = rownames(otu_matrix))
  }

  # --- Extract Sample Metadata (colData) ---
  sample_data_df <- tryCatch(
    as.data.frame(microbiome::meta(physeq)),
    error = function(e) NULL
  )

  if (!is.null(sample_data_df) && ncol(sample_data_df) > 0) {
    rownames(sample_data_df) <- phyloseq::sample_names(physeq)
  } else {
    sample_data_df <- NULL
    warning("Sample metadata is missing. Proceeding without colData.")
  }

  # --- Extract Phylogenetic Tree (if available) ---
  tree <- tryCatch(
    phyloseq::phy_tree(physeq),
    error = function(e) NULL
  )

  if (!is.null(tree) && !ape::is.rooted(tree)) {
    message("Warning: Phylogenetic tree is unrooted. Proceeding without modification.")
    tree <- ape::as.phylo(tree)  # Ensure compatibility
  }

  # --- Extract Reference Sequences (if available) ---
  ref_sequences <- tryCatch(
    phyloseq::refseq(physeq),
    error = function(e) NULL
  )

  # Ensure refseq row names match the OTU table row names
  if (!is.null(ref_sequences)) {
    refseq_df <- data.frame(RefSeq = as.character(ref_sequences), row.names = names(ref_sequences))
    refseq_df <- refseq_df[rownames(otu_matrix), , drop = FALSE]
  } else {
    refseq_df <- NULL
  }

  # --- Construct TreeSummarizedExperiment ---
  tse_args <- list(
    assays = list(counts = otu_matrix),
    rowData = tax_data
  )

  # Include sample metadata only if it exists
  if (!is.null(sample_data_df)) {
    tse_args$colData <- sample_data_df
  }

  # Include tree only if it exists
  if (!is.null(tree)) {
    tse_args$rowTree <- tree
  }

  tse <- do.call(TreeSummarizedExperiment::TreeSummarizedExperiment, tse_args)

  # Corrected: Store reference sequences using metadata from S4Vectors
  if (!is.null(refseq_df)) {
    S4Vectors::metadata(tse)$refseq <- refseq_df
  }

  return(tse)
}


# Usage Example:
# tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

# Print to confirm structure
#print(tse_16SOTU)
