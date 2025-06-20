#' @title Convert Relative ASV/OTU Counts to Absolute Counts
#'
#' @description Converts relative ASV counts in a `phyloseq` or `TreeSummarizedExperiment` (TSE) object
#' to absolute counts by multiplying ASV counts by provided scaling factors.
#' Ensures phylogenetic tree (`rowTree`) and reference sequences (`refseq`) are **retained** if possible.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @param scaling_factors A named numeric vector of scaling factors for each sample.
#' @param output_dir A character string specifying the directory to save the output file (default: `NULL`).
#' @return A list containing:
#'  - `absolute_counts`: A data frame of absolute counts.
#'  - `obj_adj`: The modified `phyloseq` or `TreeSummarizedExperiment` object.
#' @note The output file `absolute_counts.csv` is written to the specified directory. For CRAN compliance, use `tempdir()` when saving files inside examples or vignettes.
#' @importFrom phyloseq phyloseq otu_table tax_table sample_data phy_tree refseq
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom ape drop.tip
#' @importFrom utils write.csv
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   spiked_cells <- 1847
#'   species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- "Tetragenococcus_halophilus"
#'
#'   spiked_16S_OTU <- phyloseq::subset_samples(physeq_16SOTU, spiked.volume %in% c("2", "1"))
#'   Spiked_16S_sum_scaled <- Pre_processing_species(
#'     spiked_16S_OTU,
#'     species_name,
#'     merge_method = "sum",
#'     output_file = file.path(tempdir(), "merged_physeq_sum.rds")
#'   )
#'
#'   result <- calculate_spikeIn_factors(
#'     Spiked_16S_sum_scaled,
#'     spiked_cells,
#'     merged_spiked_species
#'   )
#'   scaling_factors <- result$scaling_factors
#'
#'   result_physeq <- convert_to_absolute_counts(
#'     Spiked_16S_sum_scaled,
#'     scaling_factors,
#'     output_dir = tempdir()
#'   )
#'   abs_counts_physeq <- result_physeq$absolute_counts
#'   physeq_adj <- result_physeq$obj_adj
#'
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   spiked_16S_OTU_TSE <- tse_16SOTU[, tse_16SOTU$spiked.volume %in% c("2", "1")]
#'
#'   Spiked_16S_sum_tse_scaled <- Pre_processing_species(
#'     spiked_16S_OTU_TSE,
#'     species_name,
#'     merge_method = "sum",
#'     output_file = file.path(tempdir(), "merged_tse_sum.rds")
#'   )
#'
#'   result <- calculate_spikeIn_factors(
#'     Spiked_16S_sum_tse_scaled,
#'     spiked_cells,
#'     merged_spiked_species
#'   )
#'   scaling_factors <- result$scaling_factors
#'
#'   result_tse <- convert_to_absolute_counts(
#'     spiked_16S_OTU_TSE,
#'     scaling_factors,
#'     output_dir = tempdir()
#'   )
#'   abs_counts_tse <- result_tse$absolute_counts
#'   tse_adj <- result_tse$obj_adj
#' }
#' @export
convert_to_absolute_counts <- function(obj, scaling_factors, output_dir = NULL) {
  # ==============================
  # Validate Object Type
  # ==============================
  is_physeq <- inherits(obj, "phyloseq")
  is_tse <- inherits(obj, "TreeSummarizedExperiment")

  if (!is_physeq && !is_tse) {
    stop("Input must be a `phyloseq` or `TreeSummarizedExperiment` object.")
  }

  message("Extracting OTU Table...")
  otu_table_data <- get_otu_table(obj)

  message("Extracting Taxonomy Table...")
  tax_data <- get_tax_table(obj)

  message("Extracting Sample Metadata...")
  sample_metadata <- get_sample_data(obj)

  message("Extracting Phylogenetic Tree (if available)...")
  phy_tree <- tryCatch(if (is_physeq) phyloseq::phy_tree(obj) else TreeSummarizedExperiment::rowTree(obj), error = function(e) NULL)

  message("Extracting Reference Sequences (if available)...")
  ref_sequences <- tryCatch(if (is_physeq) phyloseq::refseq(obj) else S4Vectors::metadata(obj)$refseq, error = function(e) NULL)

  # ==============================
  # Ensure Sample Names Match Scaling Factors
  # ==============================
  otu_sample_names <- colnames(otu_table_data)

  if (!all(names(scaling_factors) %in% otu_sample_names)) {
    stop("Error: Some sample names in scaling_factors do not match OTU table sample names.")
  }

  # Ensure correct order of scaling factors
  scaling_factors <- scaling_factors[otu_sample_names]

  # ==============================
  # Convert to Absolute Counts
  # ==============================
  message("Converting to Absolute Counts...")

  # Ensure scaling factors are non-zero (replace 0s with 1)
  scaling_factors[scaling_factors == 0] <- 1

  # Multiply each OTU count by the correct scaling factor
  otu_table_abs <- sweep(otu_table_data, 2, scaling_factors, "*")
  otu_table_abs <- round(otu_table_abs)

  # Ensure non-negative counts
  otu_table_abs[otu_table_abs < 0 | is.na(otu_table_abs)] <- 0

  # ==============================
  # Prune Phylogenetic Tree to Match Remaining Taxa
  # ==============================
  if (!is.null(phy_tree)) {
    common_tips <- intersect(phy_tree$tip.label, rownames(otu_table_abs))
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

  # ==============================
  # Prune Reference Sequences to Match Remaining Taxa
  # ==============================
  if (!is.null(ref_sequences)) {
    common_seqs <- intersect(names(ref_sequences), rownames(otu_table_abs))
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

  # ==============================
  # Reconstruct Object with Absolute Counts
  # ==============================
  message("Reconstructing Object...")

  if (is_physeq) {
    obj_adj <- phyloseq::phyloseq(
      phyloseq::otu_table(otu_table_abs, taxa_are_rows = TRUE),
      if (!is.null(tax_data)) phyloseq::tax_table(as.matrix(tax_data)) else NULL,
      if (!is.null(sample_metadata)) phyloseq::sample_data(sample_metadata) else NULL,
      if (!is.null(phy_tree)) phy_tree else NULL,
      if (!is.null(ref_sequences)) phyloseq::refseq(ref_sequences) else NULL
    )
  } else if (is_tse) {
    obj_adj <- TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = otu_table_abs),
      rowData = if (!is.null(tax_data)) tax_data else NULL,
      colData = if (!is.null(sample_metadata)) sample_metadata else NULL,
      rowTree = if (!is.null(phy_tree)) phy_tree else NULL,
      metadata = if (!is.null(ref_sequences)) list(refseq = ref_sequences) else list()
    )
  }

  # ==============================
  # Save Absolute Counts
  # ==============================
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }

  output_file <- file.path(output_dir, "absolute_counts.csv")
  utils::write.csv(otu_table_abs, file = output_file, row.names = TRUE)
  message("Absolute count data saved to: ", output_file)

  # ==============================
  # Return Absolute Counts & Modified Object
  # ==============================
  return(list(
    absolute_counts = as.data.frame(otu_table_abs),
    obj_adj = obj_adj
  ))
}

# Example usage with phyloseq object
# result_physeq <- convert_to_absolute_counts(physeq_16SOTU, scaling_factors)
# abs_counts_physeq <- result_physeq$absolute_counts
# physeq_adj <- result_physeq$obj_adj

# Example usage with TreeSummarizedExperiment object
# tse_obj <- convert_phyloseq_to_tse(physeq_16SOTU)
# result_tse <- convert_to_absolute_counts(tse_obj, scaling_factors)
# abs_counts_tse <- result_tse$absolute_counts
# tse_adj <- result_tse$obj_adj
# 
