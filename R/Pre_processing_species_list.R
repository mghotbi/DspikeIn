#' @title Preprocess and Merge Spike-in Species in a Phyloseq or TSE Object
#'
#' @description
#' Merges ASVs belonging to user-defined spike-in species by summing or selecting maximum counts,
#' while preserving all available metadata (taxonomy, sample data, tree, and reference sequences).
#' This function works for both `phyloseq` and `TreeSummarizedExperiment` objects.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param spiked_species Character vector of species names to be processed (match `Species` column).
#' @param merge_method Either `"sum"` (default) to sum counts across ASVs or `"max"` to keep the most abundant ASV.
#' @param output_file Optional. Character string specifying the path to save the merged object as `.rds`.
#'
#' @return A merged object of the same class as input (`phyloseq` or `TreeSummarizedExperiment`).
#'
#' @importFrom phyloseq otu_table tax_table phy_tree refseq sample_data prune_taxa
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment assay rowData assays colData
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom ape drop.tip
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("tse", package = "DspikeIn")
#'   data("physeq", package = "DspikeIn")
#'
#'   spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#'
#'   merged_TSE <- Pre_processing_species_list(
#'     tse,
#'     spiked_species = spiked_species,
#'     merge_method = "sum"
#'   )
#'
#'   merged_physeq <- Pre_processing_species_list(
#'     physeq,
#'     spiked_species = spiked_species,
#'     merge_method = "sum"
#'   )
#' }
#'
#' @export
Pre_processing_species_list <- function(obj,
                                        spiked_species,
                                        merge_method = c("sum", "max"),
                                        output_file = NULL) {
  merge_method <- match.arg(merge_method)
  message("\u25B6 Starting pre-processing...")

  is_physeq <- inherits(obj, "phyloseq")
  is_tse <- inherits(obj, "TreeSummarizedExperiment")

  if (!is_physeq && !is_tse) stop("Input object must be a 'phyloseq' or 'TreeSummarizedExperiment'.")

  # --- Accessors ---
  otu_table_data <- get_otu_table(obj)
  tax_data <- as.data.frame(get_tax_table(obj))
  sample_metadata <- get_sample_data(obj)

  # --- Optional Tree ---
  phy_tree <- tryCatch(
    if (is_physeq) phyloseq::phy_tree(obj) else TreeSummarizedExperiment::rowTree(obj),
    error = function(e) NULL
  )

  # --- Optional RefSeq ---
  ref_sequences <- tryCatch(
    if (is_physeq) phyloseq::refseq(obj) else S4Vectors::metadata(obj)$refseq,
    error = function(e) NULL
  )

  # --- Check Taxonomy ---
  if (!"Species" %in% colnames(tax_data)) stop("'Species' column not found in taxonomy table.")

  tax_data$Species <- as.character(tax_data$Species)

  # --- Process each species ---
  for (species in spiked_species) {
    message("   > Processing: ", species)

    species_asvs <- rownames(tax_data)[tax_data$Species == species]

    if (length(species_asvs) > 1) {
      message("     Merging ", length(species_asvs), " ASVs for: ", species)

      if (merge_method == "sum") {
        sum_abundances <- colSums(otu_table_data[species_asvs, , drop = FALSE])
        otu_table_data[species_asvs[1], ] <- sum_abundances
      } else {
        max_asv <- species_asvs[which.max(rowSums(otu_table_data[species_asvs, , drop = FALSE]))]
        otu_table_data[species_asvs[1], ] <- otu_table_data[max_asv, ]
      }

      # Remove redundant ASVs
      keep <- setdiff(rownames(otu_table_data), species_asvs[-1])
      otu_table_data <- otu_table_data[keep, , drop = FALSE]
      tax_data <- tax_data[keep, , drop = FALSE]
    } else if (length(species_asvs) == 1) {
      message("     Only one ASV found for: ", species)
    } else {
      warning("No ASVs found for species: ", species)
    }
  }

  # --- Prune Tree ---
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_data))
    if (length(common_tips) > 1) {
      phy_tree <- ape::drop.tip(phy_tree, setdiff(phy_tree$tip.label, common_tips))
      message("   > Tree pruned to match taxa.")
    } else {
      phy_tree <- NULL
      message("   > Tree removed due to insufficient taxa.")
    }
  }

  # --- Prune RefSeq ---
  if (!is.null(ref_sequences)) {
    common_seqs <- intersect(names(ref_sequences), rownames(otu_table_data))
    if (length(common_seqs) > 1) {
      ref_sequences <- ref_sequences[common_seqs]
      message("   > Reference sequences pruned.")
    } else {
      ref_sequences <- NULL
      message("   > Reference sequences removed due to insufficient taxa.")
    }
  }

  # --- Reconstruct ---
  if (is_physeq) {
    components <- list(
      phyloseq::otu_table(otu_table_data, taxa_are_rows = TRUE),
      phyloseq::tax_table(as.matrix(tax_data)),
      phyloseq::sample_data(sample_metadata)
    )
    if (!is.null(phy_tree)) components <- append(components, list(phy_tree))
    if (!is.null(ref_sequences)) components <- append(components, list(phyloseq::refseq(ref_sequences)))
    obj <- do.call(phyloseq::phyloseq, components)
  } else {
    obj <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = otu_table_data),
      rowData = tax_data,
      colData = sample_metadata,
      rowTree = if (!is.null(phy_tree)) phy_tree else NULL,
      metadata = if (!is.null(ref_sequences)) list(refseq = ref_sequences) else list()
    )
  }

  # --- Optional Save ---
  if (!is.null(output_file)) {
    saveRDS(obj, file = output_file)
    message("\u2713 Merged object saved to: ", output_file)
  }

  message("\u2713 Pre-processing complete.")
  return(obj)
}


# Usage Example,
# spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
# data("tse",package = "DspikeIn")
# data("physeq",package = "DspikeIn")

# merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
# merged_physeq_sum <- Pre_processing_species_list(tse, spiked_species, merge_method = "sum")
