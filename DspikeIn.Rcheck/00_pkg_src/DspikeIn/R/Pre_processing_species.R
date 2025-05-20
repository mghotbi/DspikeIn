#' @title Pre-process taxa in a phyloseq or TSE object by merging ASVs/OTUs
#'
#' @description
#' Merges ASVs/OTUs while ensuring that the phylogenetic tree and reference sequences remain intact.
#' The provided taxonomic name(s) will be searched across **all taxonomic levels** (e.g., Kingdom, Phylum, ..Genus, Species).
#' If tree or refseq become mismatched, they are pruned or removed safely.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param species_name A character vector of **exact** taxonomic names to merge (matched across all taxonomy levels).
#' @param merge_method Method used to merge counts: `"sum"` (default) or `"max"`.
#' @param output_file Optional file path to save the processed object (e.g., `file.path(tempdir(), "output.rds")`).
#' @return A processed `phyloseq` or `TreeSummarizedExperiment` object with merged ASVs/OTUs.
#'
#' @importFrom phyloseq otu_table tax_table sample_data phy_tree phyloseq refseq
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom ape drop.tip
#' @examples
#' library(DspikeIn)
#' data("physeq_16SOTU", package = "DspikeIn")
#'
#' species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'
#' # Merge species in phyloseq format
#' merged_sum <- Pre_processing_species(
#'   physeq_16SOTU,
#'   species_name,
#'   merge_method = "sum"
#' )
#'
#' # Convert phyloseq to TSE format
#' tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#' # Merge species in TSE format and write to tempdir
#' output_rds <- file.path(tempdir(), "merged_TSE_sum.rds")
#'
#' merged_TSE_sum <- Pre_processing_species(
#'   tse_16SOTU,
#'   species_name,
#'   merge_method = "sum",
#'   output_file = output_rds
#' )
#' @export
Pre_processing_species <- function(obj, species_name, merge_method = c("sum", "max"), output_file = NULL) {
  merge_method <- match.arg(merge_method)
  message("Starting pre-processing...")

  # Detect object type
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
  ref_sequences <- tryCatch(if (is_physeq) phyloseq::refseq(obj) else S4Vectors::metadata(obj)$refseq, error = function(e) NULL)

  # Convert all taxonomy columns to character
  tax_data[] <- lapply(tax_data, as.character)
  message("Checking taxonomy table...")

  # Process each species
  for (species in species_name) {
    message("Processing taxon: ", species)

    # Exact matching across ALL taxonomy levels
    species_asvs <- rownames(tax_data)[apply(tax_data, 1, function(x) any(x %in% species))]

    if (length(species_asvs) == 0) {
      warning("No ASVs found matching exactly: ", species)
      next
    }

    if (length(species_asvs) > 1) {
      message("Merging ", length(species_asvs), " ASVs for: ", species)

      if (merge_method == "sum") {
        sum_abundances <- colSums(otu_table_data[species_asvs, , drop = FALSE])
        otu_table_data[species_asvs[1], ] <- sum_abundances
      } else if (merge_method == "max") {
        max_abundance_asv <- which.max(rowSums(otu_table_data[species_asvs, , drop = FALSE]))
        otu_table_data[species_asvs[1], ] <- otu_table_data[species_asvs[max_abundance_asv], ]
      }

      # Remove redundant ASVs
      otu_table_data <- otu_table_data[setdiff(rownames(otu_table_data), species_asvs[-1]), , drop = FALSE]
      tax_data <- tax_data[setdiff(rownames(tax_data), species_asvs[-1]), , drop = FALSE]

      message("Merging completed for: ", species)
    } else {
      message("Single ASV found for: ", species, ". No merging needed.")
    }
  }

  # Prune Phylogenetic Tree to Match Remaining Taxa
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_data))
    if (length(common_tips) < length(phy_tree$tip.label)) {
      if (length(common_tips) > 1) {
        phy_tree <- ape::drop.tip(phy_tree, setdiff(phy_tree$tip.label, common_tips))
        message("Pruned phylogenetic tree to match taxa.")
      } else {
        warning("Too few taxa left in the tree after pruning. Removing tree.")
        phy_tree <- NULL
      }
    }
  }

  # Prune Reference Sequences to Match Remaining Taxa
  if (!is.null(ref_sequences)) {
    common_seqs <- intersect(names(ref_sequences), rownames(otu_table_data))
    if (length(common_seqs) < length(ref_sequences)) {
      if (length(common_seqs) > 1) {
        ref_sequences <- ref_sequences[common_seqs]
        message("Pruned reference sequences to match taxa.")
      } else {
        warning("Too few reference sequences left after pruning. Removing refseq.")
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

  # Save if requested
  if (!is.null(output_file)) {
    saveRDS(obj, file = output_file)
    message("Merged object saved to: ", output_file)
  }

  message("Pre-processing complete.")
  return(obj)
}


# Usage Example
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
# )
