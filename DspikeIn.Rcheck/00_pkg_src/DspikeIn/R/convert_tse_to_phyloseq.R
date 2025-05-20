#' @title Convert a `TreeSummarizedExperiment` to a `phyloseq` Object
#' @description Converts a `TreeSummarizedExperiment` (TSE) object back
#'              into a `phyloseq` object, ensuring the phylogenetic tree and
#'              reference sequences are restored, if present.
#'
#' @param tse A `TreeSummarizedExperiment` object containing:
#'   \itemize{
#'     \item **OTU Table (`counts`)**.
#'     \item **Taxonomy Table (`rowData`)**.
#'     \item **Sample Metadata (`colData`)**.
#'     \item **Phylogenetic Tree (`rowTree`, optional)**.
#'     \item **Reference Sequences (`metadata`, optional)**.
#'   }
#' @return A `phyloseq` object containing all preserved components.
#'
#' @importFrom phyloseq phyloseq otu_table tax_table sample_data phy_tree refseq taxa_names sample_names
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom S4Vectors metadata
#' @importFrom Biostrings DNAStringSet
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Convert phyloseq to TSE
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Convert TSE back to phyloseq
#'   phy_M <- convert_tse_to_phyloseq(tse_16SOTU)
#'   print(phy_M)
#' }
#' @export
convert_tse_to_phyloseq <- function(tse) {
  if (!inherits(tse, "TreeSummarizedExperiment")) {
    stop("Error: Input must be a valid 'TreeSummarizedExperiment' object.")
  }

  message("Extracting OTU table...")
  otu_table_matrix <- SummarizedExperiment::assay(tse, "counts")

  # Extract sample metadata
  sample_metadata <- as.data.frame(SummarizedExperiment::colData(tse))

  if (nrow(sample_metadata) == 0) {
    stop("Error: Sample metadata (colData) is empty.")
  }

  # Ensure sample names are characters
  rownames(sample_metadata) <- as.character(rownames(sample_metadata))

  # Ensure sample names in OTU table match colData row names
  if (!identical(colnames(otu_table_matrix), rownames(sample_metadata))) {
    message("Sample names mismatch detected. Aligning sample names...")
    colnames(otu_table_matrix) <- rownames(sample_metadata)
  }

  # Convert to phyloseq OTU table
  otu_table_phy <- phyloseq::otu_table(as.matrix(otu_table_matrix), taxa_are_rows = TRUE)

  # Extract and convert taxonomy table
  tax_table_df <- SummarizedExperiment::rowData(tse)
  tax_table_phy <- if (!is.null(tax_table_df) && nrow(tax_table_df) > 0) {
    phyloseq::tax_table(as.matrix(tax_table_df))
  } else {
    NULL
  }

  # Convert sample metadata to phyloseq format
  sample_data_phy <- phyloseq::sample_data(sample_metadata)

  # Ensure all sample names match before proceeding
  if (!identical(phyloseq::sample_names(otu_table_phy), phyloseq::sample_names(sample_data_phy))) {
    message("Warning: Sample names still do not match. Attempting forced alignment...")
    rownames(sample_data_phy) <- phyloseq::sample_names(otu_table_phy)
  }

  # Construct the phyloseq object
  physeq_obj <- phyloseq::phyloseq(otu_table_phy, sample_data_phy)

  if (!is.null(tax_table_phy)) {
    physeq_obj <- phyloseq::merge_phyloseq(physeq_obj, tax_table_phy)
  }

  # Extract and attach phylogenetic tree (if available)
  tree_phy <- TreeSummarizedExperiment::rowTree(tse)
  if (!is.null(tree_phy)) {
    physeq_obj <- phyloseq::merge_phyloseq(physeq_obj, phyloseq::phy_tree(tree_phy))
  }

  # Retrieve and attach reference sequences (if available)
  ref_sequences <- S4Vectors::metadata(tse)$refseq

  if (!is.null(ref_sequences)) {
    message("Checking refseq format...")

    if (is.list(ref_sequences)) {
      message("Converting refseq from list to character vector...")
      ref_sequences <- unlist(ref_sequences)
    }

    if (is.matrix(ref_sequences)) {
      message("Converting refseq from matrix to character vector...")
      ref_sequences <- as.character(ref_sequences)
    }

    if (is.character(ref_sequences)) {
      message("Converting refseq to DNAStringSet format...")
      refseq_dna <- Biostrings::DNAStringSet(ref_sequences)

      # Ensure refseq names match taxa names
      names(refseq_dna) <- taxa_names(physeq_obj)

      if (!identical(names(refseq_dna), taxa_names(physeq_obj))) {
        message("refseq names do not match taxa_names. Attempting forced alignment...")
        refseq_dna <- refseq_dna[taxa_names(physeq_obj)]
      }

      # Attach refseq to phyloseq object
      physeq_obj <- phyloseq::merge_phyloseq(physeq_obj, refseq_dna)
      message("Successfully added refseq to phyloseq object.")
    } else {
      message("Refseq format unknown. Skipping refseq assignment.")
    }
  } else {
    message("No valid reference sequences found. Proceeding without refseq.")
  }

  return(physeq_obj)
}


# Usage Example
# Convert phyloseq to TSE
# M_TSE <- convert_phyloseq_to_tse(physeq_16SOTU)

# Convert TSE back to phyloseq
# phy_M <- convert_tse_to_phyloseq(M_TSE)

# Print to confirm
# print(phy_M)
