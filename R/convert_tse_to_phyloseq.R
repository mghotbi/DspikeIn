#' @title Convert a `TreeSummarizedExperiment` to a `phyloseq` Object
#'
#' @description Converts a `TreeSummarizedExperiment` (TSE) object into a `phyloseq` object,
#' preserving key components such as the OTU table, taxonomy, sample metadata, phylogenetic tree,
#' and reference sequences (if present). This enables seamless interoperability between Bioconductor
#' and `phyloseq` workflows.
#'
#' @param tse A `TreeSummarizedExperiment` object, expected to contain:
#'   \itemize{
#'     \item OTU abundance matrix (assay named `"counts"`).
#'     \item Taxonomy data (as `rowData`).
#'     \item Sample metadata (as `colData`).
#'     \item Optional phylogenetic tree (from `rowTree()`).
#'     \item Optional reference sequences (`referenceSeq()` accessor).
#'   }
#'
#' @return A `phyloseq` object with all available components (`otu_table`, `tax_table`,
#'         `sample_data`, `phy_tree`, `refseq`) populated accordingly.
#'
#' @importFrom phyloseq phyloseq otu_table tax_table sample_data phy_tree refseq taxa_names sample_names merge_phyloseq
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom TreeSummarizedExperiment rowTree referenceSeq
#' @importFrom Biostrings DNAStringSet
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Convert to TSE and back to phyloseq
#'   tse <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   phy <- convert_tse_to_phyloseq(tse)
#'   print(phy)
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
    message("Warning: Sample names do not match. Attempting forced alignment...")
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
  ref_sequences <- tryCatch(
    TreeSummarizedExperiment::referenceSeq(tse),
    error = function(e) NULL
  )
  
  if (!is.null(ref_sequences)) {
    if (!is(ref_sequences, "DNAStringSet")) {
      message("Reference sequences found but not in DNAStringSet format. Attempting coercion...")
      ref_sequences <- Biostrings::DNAStringSet(as.character(ref_sequences))
    }
    
    # Ensure names match taxa
    if (!identical(names(ref_sequences), phyloseq::taxa_names(physeq_obj))) {
      message("Aligning refseq names to match taxa names...")
      names(ref_sequences) <- phyloseq::taxa_names(physeq_obj)
    }
    
    # Attach to phyloseq
    physeq_obj <- phyloseq::merge_phyloseq(physeq_obj, phyloseq::refseq(ref_sequences))
    message("Successfully attached reference sequences.")
  } else {
    message("No reference sequences found in TSE.")
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
# 
