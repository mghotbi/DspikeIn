#' @title Convert a `phyloseq` Object to a `TreeSummarizedExperiment`
#' 
#' @description
#' Converts a `phyloseq` object into a `TreeSummarizedExperiment` (TSE),
#' preserving key biological data components. The function supports retention of:
#' \itemize{
#'   \item OTU abundance matrix (as assay).
#'   \item Taxonomic classifications (as `rowData`).
#'   \item Sample metadata (as `colData`).
#'   \item Phylogenetic tree (as `rowTree`, if available).
#'   \item Reference sequences (as `referenceSeq`, if available).
#' }
#' This allows seamless interoperability between `phyloseq` and Bioconductor ecosystems.
#'
#' @param physeq A valid `phyloseq` object.
#'
#' @return A `TreeSummarizedExperiment` object with one or more of the following slots:
#'   \itemize{
#'     \item \code{assays}: OTU count matrix.
#'     \item \code{rowData}: Taxonomy table.
#'     \item \code{colData}: Sample metadata.
#'     \item \code{rowTree}: Phylogenetic tree (if present).
#'     \item \code{referenceSeq}: Reference sequences (if present).
#'   }
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Create a small subset for fast execution
#'   physeq_sub <- phyloseq::prune_taxa(
#'     phyloseq::taxa_names(physeq_16SOTU)[1:10],
#'     phyloseq::prune_samples(
#'       phyloseq::sample_names(physeq_16SOTU)[1:5],
#'       physeq_16SOTU
#'     )
#'   )
#'
#'   # Example transformation
#'   tse_sub <- convert_phyloseq_to_tse(physeq_sub)
#'   tse_sub
#' }
#'
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
    tree <- ape::as.phylo(tree) # Ensure compatibility
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

  if (!is.null(sample_data_df)) {
    tse_args$colData <- sample_data_df
  }
  if (!is.null(tree)) {
    tse_args$rowTree <- tree
  }
  if (!is.null(ref_sequences)) {
    tse_args$referenceSeq <- ref_sequences
  }

  tse <- do.call(TreeSummarizedExperiment::TreeSummarizedExperiment, tse_args)
  return(tse)
}


# Usage Example:
# tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

# Print to confirm structure
# print(tse_16SOTU)
