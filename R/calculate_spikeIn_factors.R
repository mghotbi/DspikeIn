#' @title Calculate Scaling Factors for Spiked Species in Phyloseq or TSE Object
#'
#' @description This function calculates scaling factors for specified spiked species
#' in a `phyloseq` or `TreeSummarizedExperiment (TSE)` object.
#' It removes spiked species, merges them, computes scaling factors,
#' applies them to the OTU table, and returns bias-corrected absolute counts along with relevant outputs.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param spiked_cells A numeric value specifying the number of spiked cells.
#' @param merged_spiked_species A character vector of spiked species to merge in the object.
#' @param output_path A character string specifying the path to save the DOCX output file.
#' Default is `NULL`, which saves the report as `"spikeIn_factors_summary.docx"`.
#'
#' @return A list containing:
#' \item{scaling_factors}{A numeric vector of calculated scaling factors, with a default factor of 1 for samples where no spiked species were detected (to prevent division by zero).}
#' \item{filtered_obj}{Filtered `phyloseq` or `TreeSummarizedExperiment` object with spiked species removed.}
#' \item{spiked_16S_total_reads}{A data frame containing total reads per sample.}
#' \item{spiked_species}{The extracted spiked species from the input object.}
#' \item{spiked_species_merged}{A merged version of the spiked species.}
#' \item{spiked_species_reads}{A data frame with read counts for merged spiked species.}
#' \item{absolute_counts}{The bias-corrected OTU table converted from relative to absolute counts.}
#' \item{normalized_obj}{The corrected `phyloseq` or `TreeSummarizedExperiment` object with all original components intact.}
#'
#' @details
#' - The metadata must contain a column named `"spiked.volume"` with appropriate values.
#' - **If no spiked species are detected in a sample, a default scaling factor of 1 is applied** (to prevent division by zero).
#' - If `TreeSummarizedExperiment` is entered and `rowTree(obj)` is invalid, the tree will be omitted.
#'
#' @importFrom phyloseq prune_taxa merge_taxa taxa_names sample_data sample_sums phy_tree otu_table tax_table phyloseq
#' @importFrom SummarizedExperiment assay rowData colData SummarizedExperiment
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom utils write.csv
#' @importFrom dplyr case_when
#' @importFrom ape drop.tip
#' @importFrom S4Vectors metadata
#' @examples
#' \donttest{
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   spiked_cells <- 1847
#'   species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- "Tetragenococcus_halophilus"
#'
#'   # Subset the phyloseq object based on spiked.volume
#'   spiked_16S_OTU <- phyloseq::subset_samples(physeq_16SOTU, spiked.volume %in% c("2", "1"))
#'
#'   # Pre-process species data by merging using the sum method
#'   merged_physeq_sum <- Pre_processing_species(
#'     spiked_16S_OTU,
#'     species_name = merged_spiked_species,
#'     merge_method = "sum",
#'     output_file = "merged_physeq_sum.rds"
#'   )
#'
#'   # Calculate spike-in factors for phyloseq object
#'   result_physeq <- calculate_spikeIn_factors(
#'   merged_physeq_sum,
#'   spiked_cells = 1874,
#'   merged_spiked_species = merged_spiked_species
#'   )
#'   scaling_factors_physeq <- result_physeq$scaling_factors
#'   print(scaling_factors_physeq)
#' }
#' @export
calculate_spikeIn_factors <- function(obj, spiked_cells, merged_spiked_species, output_path = "spikeIn_factors_output") {
  suppressMessages({

    # Ensure output directory exists
    if (!dir.exists(output_path)) {
      if (!dir.create(output_path)) stop("Failed to create directory:", output_path)
    }

    # Detect input format
    is_tse <- inherits(obj, "TreeSummarizedExperiment")
    if (is_tse) {
      cat("\U0001F500 Converting TreeSummarizedExperiment to phyloseq...\n")
      obj <- convert_tse_to_phyloseq(obj)
    }

    # Preserve original components
    tree <- tryCatch(phyloseq::phy_tree(obj), error = function(e) NULL)
    tax_data <- tryCatch(phyloseq::tax_table(obj), error = function(e) NULL)
    metadata <- tryCatch(phyloseq::sample_data(obj), error = function(e) NULL)

    cat("Extracting taxonomy and sample data...\n")

    # Ensure required columns exist
    if (!"Species" %in% colnames(tax_data)) stop("Error: 'Species' column not found in taxonomy table.")
    if (!"spiked.volume" %in% colnames(metadata)) stop("Error: 'spiked.volume' column not found in metadata.")

    cat("Removing spiked species...\n")

    # Remove spiked species
    no_spiked <- phyloseq::prune_taxa(!tax_data[, "Species"] %in% merged_spiked_species, obj)
    saveRDS(no_spiked, file.path(output_path, "Filtered_Object.rds"))

    cat("\U0001F9EE Calculating total reads per sample...\n")

    # Compute total reads
    otu_table <- phyloseq::otu_table(obj)
    total_reads <- data.frame(Sample = colnames(otu_table), Total_Reads = colSums(otu_table))
    utils::write.csv(total_reads, file.path(output_path, "Total_Reads.csv"))

    cat("\U00002796 Extracting spiked species...\n")

    # Extract and merge spiked species
    spiked_species <- phyloseq::prune_taxa(tax_data[, "Species"] %in% merged_spiked_species, obj)

    if (phyloseq::ntaxa(spiked_species) == 0) {
      warning("No spiked species found. Returning default values.")

      return(list(
        scaling_factors = setNames(rep(1, nsamples(obj)), sample_names(obj)),
        filtered_obj = obj,  # Return the original object
        spiked_total_reads = NULL,
        Total_reads = total_reads,
        spiked_species_merged = NULL,
        spiked_species_reads = NULL
      ))
    }

    saveRDS(spiked_species, file.path(output_path, "Spiked_Species.rds"))

    cat("\U00002795 Merging spiked species...\n")
    spiked_species_merged <- phyloseq::merge_taxa(spiked_species, phyloseq::taxa_names(spiked_species))
    saveRDS(spiked_species_merged, file.path(output_path, "Spiked_Species_Merged.rds"))

    cat("\U0001F9EE Calculating scaling factors...\n")

    # Compute scaling factors
    spiked_otu_table <- phyloseq::otu_table(spiked_species_merged)
    spiked_species_reads <- data.frame(Sample = colnames(spiked_otu_table), Total_Reads = colSums(spiked_otu_table))
    utils::write.csv(spiked_species_reads, file.path(output_path, "Spiked_Species_Reads.csv"))

    scaling_factors <- setNames(rep(1, nrow(total_reads)), total_reads$Sample)
    nonzero_samples <- spiked_species_reads$Total_Reads > 0
    if (sum(nonzero_samples) > 0) {
      scaling_factors[nonzero_samples] <- spiked_cells / spiked_species_reads$Total_Reads[nonzero_samples] *
        dplyr::case_when(
          metadata$spiked.volume[nonzero_samples] == 0 ~ 1,
          metadata$spiked.volume[nonzero_samples] == 0.5 ~ 1/4,
          metadata$spiked.volume[nonzero_samples] == 1   ~ 1/2,
          metadata$spiked.volume[nonzero_samples] == 2   ~ 1,
          metadata$spiked.volume[nonzero_samples] == 3   ~ 1.5,
          metadata$spiked.volume[nonzero_samples] == 4   ~ 2,
          TRUE ~ 1
        )
    }

    utils::write.csv(data.frame(Sample = total_reads$Sample, Scaling_Factor = scaling_factors),
                     file.path(output_path, "Scaling_Factors.csv"))

    return(list(
      scaling_factors = scaling_factors,
      filtered_obj = no_spiked,
      spiked_total_reads = spiked_species_reads,
      Total_reads = total_reads,
      spiked_species_merged = spiked_species_merged,
      spiked_species_reads = spiked_species_reads
    ))
  })
}

# Example usage:
# # Define the spiked species
# merged_spiked_species <- c("Tetragenococcus_halophilus")
#
# # Calculate scaling factors and generate the report
# result <- calculate_spikeIn_factors(merged_physeq_sum, 1874, merged_spiked_species)
#
# # Access the results
# scaling_factors <- result$scaling_factors
#
# Define spiked species for TSE format
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# print(merged_spiked_species %in% unique(get_tax_table(merged_sum)$Species))

# Run function for a phyloseq object
# result_physeq <- calculate_spikeIn_factors(merged_TSE_sum,
# 1874, merged_spiked_species)
# scaling_factors<-result_physeq$scaling_factors


