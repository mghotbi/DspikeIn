#' @title Pre-process species in a phyloseq or TSE object by merging ASVs/OTUs
#' @description Merges ASVs/OTUs while ensuring the phylogenetic tree and reference sequences remain intact.
#' If tree or refseq become mismatched, they are pruned or removed safely.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param species_name A character vector of species names to merge.
#' @param merge_method `"sum"` (default) or `"max"`: method for merging ASV counts.
#' @param output_file Optional file path to save the processed object.
#' @return A processed `phyloseq` or `TreeSummarizedExperiment` object.
#' @importFrom phyloseq otu_table tax_table sample_data phy_tree phyloseq refseq
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom ape drop.tip
#' @examples
#' \dontrun{
#' library(DspikeIn)
#' data("physeq_16SOTU", package = "DspikeIn")
#'  spiked_cells <- 1847
#'  species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'
#'  merged_sum <- Pre_processing_species(physeq_16SOTU, species_name, merge_method = "sum")
#'
#' # Convert phyloseq to TSE format
#' tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#' species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp." )
#'
#' # Merge species in TSE format
#' merged_TSE_sum <- Pre_processing_species(
#'     tse_16SOTU,
#'     species_name,
#'     merge_method = "sum",
#'     output_file = "merged_TSE_sum.rds"
#' )
#' }
#' @export
Pre_processing_species <- function(obj, species_name, merge_method = c("sum", "max"), output_file = NULL) {
  merge_method <- match.arg(merge_method)
  message("\U0001F504 Starting pre-processing...")

  # Detect object type
  is_physeq <- inherits(obj, "phyloseq")
  is_tse <- inherits(obj, "TreeSummarizedExperiment")

  if (!is_physeq && !is_tse) {
    stop("\U0000274C Input must be `phyloseq` or `TreeSummarizedExperiment`.")
  }

  # Extract components
  otu_table_data <- get_otu_table(obj)
  tax_data <- as.data.frame(get_tax_table(obj))
  sample_metadata <- get_sample_data(obj)
  phy_tree <- tryCatch(if (is_physeq) phyloseq::phy_tree(obj) else TreeSummarizedExperiment::rowTree(obj), error = function(e) NULL)
  ref_sequences <- tryCatch(if (is_physeq) phyloseq::refseq(obj) else S4Vectors::metadata(obj)$refseq, error = function(e) NULL)

  # Convert species column to character
  tax_data$Species <- as.character(tax_data$Species)
  message("\U0001F50D Checking taxonomy table...")

  # Process each species
  for (species in species_name) {
    message("\U0001F504 Processing species: ", species)

    species_asvs <- rownames(tax_data)[which(tax_data$Species == species)]
    if (length(species_asvs) == 0) {
      species_asvs <- rownames(tax_data)[grep(species, tax_data$Species, ignore.case = TRUE)]
      if (length(species_asvs) > 0) {
        message("\U000026A0 Using partial match for species: ", species)
      }
    }

    if (length(species_asvs) > 1) {
      if (merge_method == "sum") {
        message("\U00002795 Merging ASVs by summing abundances for: ", species)
        sum_abundances <- colSums(otu_table_data[species_asvs, , drop = FALSE])
        otu_table_data[species_asvs[1], ] <- sum_abundances
      } else if (merge_method == "max") {
        message("\U00002795 Merging ASVs using 'max' method for: ", species)
        max_abundance_asv <- which.max(rowSums(otu_table_data[species_asvs, , drop = FALSE]))
        otu_table_data[species_asvs[1], ] <- otu_table_data[max_abundance_asv, ]
      }

      # Remove redundant ASVs
      otu_table_data <- otu_table_data[setdiff(rownames(otu_table_data), species_asvs[-1]), , drop = FALSE]
      tax_data <- tax_data[setdiff(rownames(tax_data), species_asvs[-1]), , drop = FALSE]

      tax_data[species_asvs[1], "Genus"] <- tax_data[species_asvs[1], "Genus"]
      tax_data[species_asvs[1], "Species"] <- species
      message("\U00002705 Merging completed for: ", species)
    }
  }

  # **Prune Phylogenetic Tree to Match Remaining Taxa**
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_data))
    if (length(common_tips) < length(phy_tree$tip.label)) {
      if (length(common_tips) > 1) {
        phy_tree <- ape::drop.tip(phy_tree, setdiff(phy_tree$tip.label, common_tips))
        message("\U00002705 Pruned phylogenetic tree to match taxa.")
      } else {
        warning("\U000026A0 Too few taxa left in the tree after pruning. Removing tree.")
        phy_tree <- NULL
      }
    }
  }

  # **Prune Reference Sequences to Match Remaining Taxa**
  if (!is.null(ref_sequences)) {
    common_seqs <- intersect(names(ref_sequences), rownames(otu_table_data))
    if (length(common_seqs) < length(ref_sequences)) {
      if (length(common_seqs) > 1) {
        ref_sequences <- ref_sequences[common_seqs]
        message("\U00002705 Pruned reference sequences to match taxa.")
      } else {
        warning("\U000026A0 Too few reference sequences left after pruning. Removing refseq.")
        ref_sequences <- NULL
      }
    }
  }

  # Reconstruct Object
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
      rowTree = if (!is.null(phy_tree)) phy_tree else NULL,
      metadata = if (!is.null(ref_sequences)) list(refseq = ref_sequences) else list()
    )
  }

  if (!is.null(output_file)) {
    saveRDS(obj, file = output_file)
    message("\U0001F4BE Merged object saved to: ", output_file)
  }

  message("\U00002705 Pre-processing complete.")
  return(obj)
}

#Usage Example
# data("physeq_16SOTU", package="DspikeIn")
# species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
# merged_sum <- Pre_processing_species(physeq_16SOTU, species_name, merge_method = "sum")

# tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
# species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
# Merge species in TSE format
# merged_TSE_sum <- Pre_processing_species(
#  tse_16SOTU,
#  species_name,
#  merge_method = "sum",
#  output_file = "merged_TSE_sum.rds"
#)




