#' @title Pre-process taxa in a phyloseq or TSE object by merging ASVs/OTUs using hashcodes
#'
#' @description
#' Merges ASVs/OTUs using a list of exact hashcodes (typically OTU IDs from QIIME2).
#' Supports both `phyloseq` and `TreeSummarizedExperiment` (TSE) objects. Ensures that
#' phylogenetic trees and reference sequences are retained and pruned as needed to remain
#' consistent with merged taxa.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param hashcodes A character vector of OTU hashcodes to be merged.
#' @param merge_method Merge strategy to apply to counts: `"sum"` (default) or `"max"`.
#' @param output_file Optional path to save the resulting object as an `.rds` file.
#'
#' @return A processed `phyloseq` or `TreeSummarizedExperiment` object with merged ASVs/OTUs.
#'
#' @importFrom phyloseq phyloseq otu_table tax_table sample_data phy_tree refseq merge_phyloseq
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree referenceSeq referenceSeq<-
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom ape drop.tip
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Example 1: phyloseq input
#'   subset_physeq <- phyloseq::subset_taxa(
#'     physeq_16SOTU,
#'     Species %in% c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   )
#'   hashcodes_phy <- rownames(phyloseq::otu_table(subset_physeq))
#'
#'   merged_physeq <- Pre_processing_hashcodes(
#'     obj = physeq_16SOTU,
#'     hashcodes = hashcodes_phy,
#'     merge_method = "sum"
#'   )
#'
#'   # Example 2: TreeSummarizedExperiment (TSE) input
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   row_annots <- SummarizedExperiment::rowData(tse_16SOTU)
#'   subset_tse <- tse_16SOTU[
#'     which(row_annots$Species %in% c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")),
#'   ]
#'   hashcodes_tse <- rownames(SummarizedExperiment::assay(subset_tse, "counts"))
#'
#'   merged_tse <- Pre_processing_hashcodes(
#'     obj = tse_16SOTU,
#'     hashcodes = hashcodes_tse,
#'     merge_method = "max"
#'   )
#'
#'   # Optional cleanup
#'   output_files <- c(
#'     file.path(tempdir(), "merged_physeq_sum_processed.rds"),
#'     file.path(tempdir(), "merged_tse_max_processed.rds")
#'   )
#'   for (f in output_files) {
#'     if (file.exists(f)) file.remove(f)
#'   }
#'
#'   rm(subset_physeq, hashcodes_phy, merged_physeq,
#'      tse_16SOTU, subset_tse, hashcodes_tse, merged_tse, output_files, row_annots)
#' }
#' @export
Pre_processing_hashcodes <- function(obj, hashcodes, merge_method = c("sum", "max"), output_file = NULL) {
  merge_method <- match.arg(merge_method)
  message("Starting pre-processing...")
  
  is_physeq <- inherits(obj, "phyloseq")
  is_tse <- inherits(obj, "TreeSummarizedExperiment")
  
  if (!is_physeq && !is_tse) {
    stop("Input must be `phyloseq` or `TreeSummarizedExperiment`.")
  }
  
  # Extract components
  otu_table_data <- get_otu_table(obj)
  tax_data <- as.data.frame(get_tax_table(obj))
  sample_metadata <- get_sample_data(obj)
  phy_tree <- tryCatch(if (is_physeq) phyloseq::phy_tree(obj) else TreeSummarizedExperiment::rowTree(obj), error = function(e) NULL)
  ref_sequences <- tryCatch(
    if (is_physeq) phyloseq::refseq(obj) else TreeSummarizedExperiment::referenceSeq(obj),
    error = function(e) NULL
  )
  
  # Validate input hashcodes
  hashcodes <- intersect(hashcodes, rownames(otu_table_data))
  if (length(hashcodes) < 2) {
    stop("At least two matching hashcodes are required for merging.")
  }
  
  # Merge selected ASVs/OTUs
  if (merge_method == "sum") {
    merged_counts <- colSums(otu_table_data[hashcodes, , drop = FALSE])
    otu_table_data[hashcodes[1], ] <- merged_counts
  } else {
    max_idx <- which.max(rowSums(otu_table_data[hashcodes, , drop = FALSE]))
    otu_table_data[hashcodes[1], ] <- otu_table_data[hashcodes[max_idx], ]
  }
  
  # Drop redundant rows
  otu_table_data <- otu_table_data[setdiff(rownames(otu_table_data), hashcodes[-1]), , drop = FALSE]
  tax_data <- tax_data[setdiff(rownames(tax_data), hashcodes[-1]), , drop = FALSE]
  
  # Prune tree if needed
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_data))
    if (length(common_tips) > 1) {
      phy_tree <- ape::drop.tip(phy_tree, setdiff(phy_tree$tip.label, common_tips))
    } else {
      warning("Too few taxa remain for tree. Tree removed.")
      phy_tree <- NULL
    }
  }
  
  # Prune refseq if needed
  if (!is.null(ref_sequences)) {
    common_names <- intersect(names(ref_sequences), rownames(otu_table_data))
    if (length(common_names) > 1) {
      ref_sequences <- ref_sequences[common_names]
    } else {
      warning("Too few sequences left. RefSeq removed.")
      ref_sequences <- NULL
    }
  }
  
  # Rebuild object
  if (is_physeq) {
    obj <- phyloseq::phyloseq(
      phyloseq::otu_table(otu_table_data, taxa_are_rows = TRUE),
      phyloseq::tax_table(as.matrix(tax_data)),
      phyloseq::sample_data(sample_metadata),
      if (!is.null(phy_tree)) phy_tree else NULL,
      if (!is.null(ref_sequences)) phyloseq::refseq(ref_sequences) else NULL
    )
  } else {
    obj <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = otu_table_data),
      rowData = tax_data,
      colData = sample_metadata,
      rowTree = if (!is.null(phy_tree)) phy_tree else NULL
    )
    
    if (!is.null(ref_sequences)) {
      TreeSummarizedExperiment::referenceSeq(obj) <- ref_sequences
    }
  }
  
  # Save if needed
  if (!is.null(output_file)) {
    saveRDS(obj, file = output_file)
    message("Saved processed object to: ", output_file)
  }
  
  message("Pre-processing complete.")
  return(obj)
}


# Example usage:
# Tetragenococcus <- phyloseq::subset_taxa(physeq_16SOTU,
# Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp.")
# hashcodes <- row.names(phyloseq::otu_table(Tetragenococcus))
# processed_data_sum <- Pre_processing_hashcodes(physeq_16SOTU, hashcodes, merge_method = "sum",
# output_prefix = "merged_physeq_sum")

# Tetragenococcus_TSE <- convert_phyloseq_to_tse(Tetragenococcus)
# hashcodes <- rownames(get_otu_table(Tetragenococcus_TSE))
# tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
# processed_data_max <- Pre_processing_hashcodes(tse_16SOTU, hashcodes, merge_method = "max",
# output_prefix = "merged_physeq_max")
