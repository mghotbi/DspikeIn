#' @title Filter and Split Abundance Data (Saves Phyloseq & TSE Separately)
#' @description Filters low-abundance taxa from a `phyloseq` or `TreeSummarizedExperiment` object,
#'              splits data into high- and low-abundance groups, and saves all results **in both formats**.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param threshold A numeric value indicating the mean abundance threshold for filtering and splitting.
#' @param output_prefix A character string specifying the filename prefix for saved results.
#'
#' @return A list containing:
#' - `$filtered`: Filtered objects in both formats.
#' - `$high`: High-abundance taxa in both formats.
#' - `$low`: Low-abundance taxa in both formats.
#'
#' @details
#' The function:
#' - Removes taxa with mean abundance below `threshold`
#' - Splits remaining data into **high-abundance** and **low-abundance** groups
#' - Saves all results as **both `phyloseq` and `TreeSummarizedExperiment` formats**
#' - Notifies users where files are saved via `cat()`
#'
#' @examples
#' \donttest{
#' data("physeq_ITSOTU", package="DspikeIn")
#' results <- filter_and_split_abundance(physeq_ITSOTU,
#' threshold = 0.05,
#' output_prefix = "abundance_analysis")
#' tse_ITSOTU <-convert_phyloseq_to_tse(physeq_ITSOTU)
#' results <- filter_and_split_abundance(tse_ITSOTU,
#' threshold = 0.05, output_prefix = "abundance_analysis")
#' }
#'
#' @importFrom phyloseq prune_taxa
#' @importFrom SummarizedExperiment assay
#' @export
filter_and_split_abundance <- function(obj, threshold = 0.01, output_prefix = "abundance_analysis") {
  # Ensure object is valid
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  # Extract OTU table
  otu <- get_otu_table(obj)
  high_otu <- otu[rowMeans(otu) > threshold, , drop = FALSE]
  low_otu <- otu[rowMeans(otu) <= threshold, , drop = FALSE]

  # Process based on object type
  if (inherits(obj, "phyloseq")) {
    filtered_obj_phy <- phyloseq::prune_taxa(rownames(high_otu), obj)
    high_obj_phy <- phyloseq::prune_taxa(rownames(high_otu), obj)
    low_obj_phy <- phyloseq::prune_taxa(rownames(low_otu), obj)

    # Convert phyloseq to TSE
    filtered_obj_tse <- convert_phyloseq_to_tse(filtered_obj_phy)
    high_obj_tse <- convert_phyloseq_to_tse(high_obj_phy)
    low_obj_tse <- convert_phyloseq_to_tse(low_obj_phy)

  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    filtered_obj_tse <- obj[rownames(high_otu), ]
    high_obj_tse <- obj[rownames(high_otu), ]
    low_obj_tse <- obj[rownames(low_otu), ]

    # Convert TSE to phyloseq
    filtered_obj_phy <- convert_tse_to_phyloseq(filtered_obj_tse)
    high_obj_phy <- convert_tse_to_phyloseq(high_obj_tse)
    low_obj_phy <- convert_tse_to_phyloseq(low_obj_tse)
  }

  # Save results with distinct filenames
  saveRDS(filtered_obj_phy, paste0(output_prefix, "_filtered_phyloseq.rds"))
  saveRDS(filtered_obj_tse, paste0(output_prefix, "_filtered_tse.rds"))
  saveRDS(high_obj_phy, paste0(output_prefix, "_high_abundance_phyloseq.rds"))
  saveRDS(high_obj_tse, paste0(output_prefix, "_high_abundance_tse.rds"))
  saveRDS(low_obj_phy, paste0(output_prefix, "_low_abundance_phyloseq.rds"))
  saveRDS(low_obj_tse, paste0(output_prefix, "_low_abundance_tse.rds"))

  # Inform user where files are saved
  cat("\n **Filtered and split datasets saved!**\n")
  cat("\U0001F5C2 Filtered (phyloseq): ", paste0(output_prefix, "_filtered_phyloseq.rds"), "\n")
  cat("\U0001F5C2 Filtered (TSE): ", paste0(output_prefix, "_filtered_tse.rds"), "\n")
  cat("\U0001F5C2 High-abundance (phyloseq): ", paste0(output_prefix, "_high_abundance_phyloseq.rds"), "\n")
  cat("\U0001F5C2 High-abundance (TSE): ", paste0(output_prefix, "_high_abundance_tse.rds"), "\n")
  cat("\U0001F5C2 Low-abundance (phyloseq): ", paste0(output_prefix, "_low_abundance_phyloseq.rds"), "\n")
  cat("\U0001F5C2 Low-abundance (TSE): ", paste0(output_prefix, "_low_abundance_tse.rds"), "\n")

  # Return results as a named list
  return(list(
    filtered = list(phyloseq = filtered_obj_phy, tse = filtered_obj_tse),
    high = list(phyloseq = high_obj_phy, tse = high_obj_tse),
    low = list(phyloseq = low_obj_phy, tse = low_obj_tse)
  ))
}

#Usage Example:
# results <- filter_and_split_abundance(physeq_ITSOTU,
# threshold = 0.05, output_prefix = "abundance_analysis")

# tse_ITSOTU <-convert_phyloseq_to_tse(physeq_ITSOTU)
# results <- filter_and_split_abundance(tse_ITSOTU,
# threshold = 0.05, output_prefix = "abundance_analysis")
