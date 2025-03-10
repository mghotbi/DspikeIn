#' @title Pre-process a list of spiked-in species in a phyloseq or TSE object
#' @description Merges ASVs based on a specified method while preserving all metadata
#' (taxonomy, sample data, phylogenetic tree, and reference sequences, if available).
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param spiked_species A character vector of species names to be processed (Genus Species format).
#' @param merge_method Either `"sum"` (sum counts) or `"max"` (keep max abundance ASV). Default is `"sum"`.
#' @param output_file Optional: File path to save the merged object as an `.rds` file.
#' @return A `phyloseq` or `TreeSummarizedExperiment` object with merged species.
#'
#' @importFrom phyloseq otu_table tax_table phy_tree refseq sample_data prune_taxa
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment assay rowData assays colData
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom ape is.rooted drop.tip
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' library(TreeSummarizedExperiment)
#'
#' data("tse", package = "DspikeIn")
#' data("physeq", package = "DspikeIn")
#'
#' # Define spiked species for merging
#' spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
#'
#' # Process the dataset and merge ASVs
#' merged_TSE <- Pre_processing_species_list(tse, spiked_species, merge_method = "sum")
#' merged_physeq <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
#' }
#' @export
Pre_processing_species_list <- function(obj, spiked_species, merge_method = c("sum", "max"), output_file = NULL) {

  merge_method <- match.arg(merge_method)
  message("Starting pre-processing for spiked-in species...")

  # Step 1: Detect input format
  is_physeq <- inherits(obj, "phyloseq")
  is_tse <- inherits(obj, "TreeSummarizedExperiment")

  if (!is_physeq && !is_tse) {
    stop("Input object must be either a `phyloseq` or `TreeSummarizedExperiment`.")
  }

  # Step 2: Retrieve components
  otu_table_data <- get_otu_table(obj)
  tax_data <- as.data.frame(get_tax_table(obj))
  sample_metadata <- get_sample_data(obj)

  # Retrieve phylogenetic tree (if available)
  phy_tree <- tryCatch(
    if (is_physeq) phyloseq::phy_tree(obj) else TreeSummarizedExperiment::rowTree(obj),
    error = function(e) NULL
  )

  # Retrieve reference sequences (if available)
  ref_sequences <- tryCatch(
    if (is_physeq) phyloseq::refseq(obj) else S4Vectors::metadata(obj)$refseq,
    error = function(e) NULL
  )

  # Step 3: Ensure `Species` column exists
  if (!"Species" %in% colnames(tax_data)) stop("Error: 'Species' column not found in taxonomy table.")
  tax_data$Species <- as.character(tax_data$Species)

  message("Checking taxonomy table...")

  # Step 4: Process each species
  for (species in spiked_species) {
    message("Processing species: ", species)

    species_asvs <- rownames(tax_data)[which(tax_data$Species == species)]

    if (length(species_asvs) > 1) {
      message("Merging ", length(species_asvs), " ASVs for species: ", species)

      if (merge_method == "sum") {
        sum_abundances <- colSums(otu_table_data[species_asvs, , drop = FALSE])
        otu_table_data[species_asvs[1], ] <- sum_abundances
      } else if (merge_method == "max") {
        max_abundance_asv <- species_asvs[which.max(rowSums(otu_table_data[species_asvs, , drop = FALSE]))]
        otu_table_data[species_asvs[1], ] <- otu_table_data[max_abundance_asv, ]
      }

      # Remove extra ASVs
      otu_table_data <- otu_table_data[setdiff(rownames(otu_table_data), species_asvs[-1]), , drop = FALSE]
      tax_data <- tax_data[setdiff(rownames(tax_data), species_asvs[-1]), , drop = FALSE]
    } else if (length(species_asvs) == 1) {
      message("No merging needed; only one ASV/OTU for species: ", species)
    } else {
      warning("Warning: No ASVs/OTUs found for species: ", species)
    }
  }

  # Step 5: Prune Tree and RefSeq for TSE
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_data))
    if (length(common_tips) < length(phy_tree$tip.label)) {
      if (length(common_tips) > 1) {
        phy_tree <- ape::drop.tip(phy_tree, setdiff(phy_tree$tip.label, common_tips))
        message("Pruned tree to match remaining taxa.")
      } else {
        warning("Tree has too few taxa after pruning. Removing tree.")
        phy_tree <- NULL
      }
    }
  }

  if (!is.null(ref_sequences)) {
    common_seqs <- intersect(names(ref_sequences), rownames(otu_table_data))
    if (length(common_seqs) < length(ref_sequences)) {
      if (length(common_seqs) > 1) {
        ref_sequences <- ref_sequences[common_seqs]
        message("Pruned reference sequences to match remaining taxa.")
      } else {
        warning("Reference sequences have too few taxa after pruning. Removing refseq.")
        ref_sequences <- NULL
      }
    }
  }

  # Step 6: Reconstruct Object
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

  if (!is.null(output_file)) {
    saveRDS(obj, file = output_file)
    message("Merged object saved to: ", output_file)
  }

  message("Pre-processing complete.")
  return(obj)
}


# Usage Example,
# spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
# data("tse",package = "DspikeIn")
# data("physeq",package = "DspikeIn")

# merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")
# merged_physeq_sum <- Pre_processing_species_list(tse, spiked_species, merge_method = "sum")

