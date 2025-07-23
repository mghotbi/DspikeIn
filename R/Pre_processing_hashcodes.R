#' @title Pre-process phyloseq or TSE object based on hashcodes
#' @description Subsets, merges, and saves taxa based on hashcodes and a specified merge method ("sum" or "max").
#' This function pre-processes a `phyloseq` or `TreeSummarizedExperiment` (TSE) object by subsetting, merging,
#' and saving taxa based on provided hashcodes. It retains taxonomic information and creates intermediate datasets
#' for further downstream analysis.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param hashcodes A character vector of taxon hashcodes (OTU row names).
#' @param merge_method The method to merge taxa: "sum" or "max".
#' @param output_prefix A prefix for the output file names.
#' @return A processed phyloseq or TSE object.
#'
#' @importFrom phyloseq taxa_names tax_table otu_table merge_taxa prune_taxa merge_phyloseq sample_data
#' @importFrom TreeSummarizedExperiment rowTree
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Subset to Tetragenococcus species
#'   tetragenococcus_physeq <- phyloseq::subset_taxa(
#'     physeq_16SOTU,
#'     Species %in% c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   )
#'
#'   # Extract OTU IDs (hashcodes) for phyloseq object
#'   hashcodes_physeq <- rownames(phyloseq::otu_table(tetragenococcus_physeq))
#'
#'   # Remove previous output file if exists
#'   if (file.exists("merged_physeq_processed.rds"))
#'     file.remove("merged_physeq_processed.rds")
#'
#'   # Run merging with "sum" method for phyloseq
#'   processed_sum <- Pre_processing_hashcodes(
#'     physeq_16SOTU,
#'     hashcodes = hashcodes_physeq,
#'     merge_method = "sum"
#'   )
#'
#'   # Convert to TreeSummarizedExperiment (TSE)
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   tetragenococcus_TSE <- convert_phyloseq_to_tse(tetragenococcus_physeq)
#'
#'   # Extract hashcodes for TSE
#'   hashcodes_tse <- rownames(tetragenococcus_TSE)
#'
#'   # Run merging with "max" method for TSE
#'   processed_max <- Pre_processing_hashcodes(
#'     tse_16SOTU,
#'     hashcodes = hashcodes_tse,
#'     merge_method = "max"
#'   )

#'
#'   # Final cleanup of written file
#'   file.remove("merged_physeq_processed.rds")
#' }
#' @export
Pre_processing_hashcodes <- function(obj, hashcodes, merge_method = c("sum", "max"), output_prefix = "merged_physeq") {
  merge_method <- match.arg(merge_method)
  message("Starting pre-processing...")

  #  Convert TSE to phyloseq if necessary
  if (inherits(obj, "TreeSummarizedExperiment")) {
    message("Converting TSE to phyloseq...")
    obj <- convert_tse_to_phyloseq(obj)
  }

  #  Extract data using accessors
  otu_matrix <- get_otu_table(obj)
  tax_table_df <- get_tax_table(obj)
  sample_metadata <- get_sample_data(obj)

  # Ensure valid data
  if (is.null(otu_matrix) || nrow(otu_matrix) == 0) {
    stop("Error: OTU table is empty or missing.")
  }

  if (!all(hashcodes %in% rownames(otu_matrix))) {
    stop("Error: One or more hashcodes not found in the dataset.")
  }

  #  Merge method: "sum" (sum abundances)
  if (merge_method == "sum") {
    message("Merging taxa using 'sum' method...")
    obj_merged <- phyloseq::merge_taxa(obj, hashcodes)
  }

  #  Merge method: "max" (keep max abundance per sample)
  else if (merge_method == "max") {
    message("Merging taxa using 'max' method...")

    # Extract and process OTU table
    otu_selected <- otu_matrix[hashcodes, , drop = FALSE]
    max_abundances <- apply(otu_selected, 2, max)
    max_hashcodes <- hashcodes[apply(otu_selected, 2, which.max)]

    # Create a new OTU table with max abundances
    new_otu_table <- otu_selected[max_hashcodes[1], , drop = FALSE]
    new_otu_table[] <- max_abundances

    # Retain taxonomic info
    new_tax_table <- tax_table_df[max_hashcodes[1], , drop = FALSE]

    #  Build a new phyloseq object with max-abundance taxa
    max_phyloseq <- phyloseq::phyloseq(
      phyloseq::otu_table(new_otu_table, taxa_are_rows = TRUE),
      phyloseq::tax_table(as.matrix(new_tax_table)),
      phyloseq::sample_data(sample_metadata)
    )

    # Remove original hashcodes from dataset
    obj_pruned <- phyloseq::prune_taxa(!rownames(otu_matrix) %in% hashcodes, obj)

    # Merge pruned object with max-abundance taxa
    obj_merged <- phyloseq::merge_phyloseq(obj_pruned, max_phyloseq)
  }

  #  Save processed data
  processed_data_path <- paste0(output_prefix, "_processed.rds")
  saveRDS(obj_merged, processed_data_path)
  message("Saved processed data file: ", processed_data_path)

  message("Pre-processing complete.")
  return(obj_merged)
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
