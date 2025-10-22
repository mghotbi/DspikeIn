#' @title Calculate Scaling Factors for Spiked Species in Phyloseq or TSE Object
#'
#' @description Calculates scaling factors for specified spike-in species or genera in
#' a `phyloseq` or `TreeSummarizedExperiment (TSE)` object. It supports genus/species-level detection,
#' removes spike-ins, merges them, computes scaling factors, and returns a bias-corrected absolute count matrix.
<<<<<<< HEAD
=======
#' The function automatically handles:Species or genus-level spike-in identification, Safe tree and taxonomy synchronization,
#' Volume-based scaling (via the spiked.volume field in metadata), 
#' Optional export of intermediate results for traceability (Total_Reads.csv, Spiked_Reads.csv, Scaling_Factors.csv).
>>>>>>> bioc-rev1
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param spiked_cells A numeric value for the number of spiked cells per unit volume.
#' @param merged_spiked_species A character vector of spiked taxon names (species or genus).
#' @param output_path Optional directory path to save intermediate files (default is NULL).
#'
#' @return A list with:
#' \item{scaling_factors}{Named numeric vector of scaling factors per sample.}
#' \item{filtered_obj}{Input object with spike-in taxa removed.}
#' \item{spiked_species_reads}{Data frame with spike-in reads per sample.}
#' \item{total_reads}{Data frame with total reads per sample.}
#' \item{spiked_species_merged}{The merged spike-in taxa as a phyloseq object.}
#' \item{tree}{Original phylogenetic tree, if available.}
#'
#' @importFrom phyloseq prune_taxa merge_taxa taxa_names sample_data sample_names sample_sums phy_tree otu_table tax_table phyloseq
#' @importFrom utils write.csv
#' @importFrom dplyr case_when
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE) &&
#'   requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   spiked_cells <- 1847
#'   species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
#'   merged_spiked_species <- "Tetragenococcus_halophilus"
#'
#'   # --- Phyloseq Example ---
#'   spiked_16S_OTU <- phyloseq::subset_samples(physeq_16SOTU, spiked.volume %in% c("2", "1"))
#'   temp_output_file <- file.path(tempdir(), "merged_physeq_sum.rds")
#'   output_dir <- file.path(tempdir(), "spikeIn_factors_output")
#'
#'   merged_physeq_sum <- Pre_processing_species(
#'     spiked_16S_OTU,
#'     species_name = merged_spiked_species,
#'     merge_method = "sum",
#'     output_file = temp_output_file
#'   )
#'
#'   result_physeq <- calculate_spikeIn_factors(
#'     merged_physeq_sum,
#'     spiked_cells = spiked_cells,
#'     merged_spiked_species = merged_spiked_species,
#'     output_path = output_dir
#'   )
#'
#'   print(result_physeq$scaling_factors)
#'
#'   if (file.exists(temp_output_file)) unlink(temp_output_file, force = TRUE)
#'   if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)
#'
#'   # --- TSE Example ---
#'   tse_data <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   merged_tse_sum <- Pre_processing_species(
#'     tse_data,
#'     species_name = merged_spiked_species,
#'     merge_method = "sum"
#'   )
#'
#'   result_tse <- calculate_spikeIn_factors(
#'     merged_tse_sum,
#'     spiked_cells = spiked_cells,
#'     merged_spiked_species = merged_spiked_species
#'   )
#'
#'   print(result_tse$scaling_factors)
#'
#'   # --- Final cleanup of any extra leftover RDS files ---
#'   leftover_rds <- list.files(tempdir(), pattern = "merged_physeq.*\\.rds$", full.names = TRUE)
#'   file.remove(leftover_rds[file.exists(leftover_rds)])
#' }
#' @export
calculate_spikeIn_factors <- function(obj, spiked_cells, merged_spiked_species, output_path = NULL) {
  if (!is.null(output_path) && !dir.exists(output_path)) {
    if (!dir.create(output_path, recursive = TRUE)) stop("Failed to create directory: ", output_path)
  }

  # Check if it's a TreeSummarizedExperiment
  is_tse <- inherits(obj, "TreeSummarizedExperiment")
  if (is_tse) {
    message("Converting TreeSummarizedExperiment to phyloseq...")
    obj <- convert_tse_to_phyloseq(obj)
  }

  # Safely extract components
  tree <- tryCatch(
    {
      phyloseq::phy_tree(obj)
    },
    error = function(e) {
      message("No phylogenetic tree found. Continuing without it.")
      NULL
    }
  )

  refseq <- tryCatch(
    {
      phyloseq::refseq(obj)
    },
    error = function(e) {
      message("No reference sequences found. Continuing without them.")
      NULL
    }
  )

  tax_data <- tryCatch(phyloseq::tax_table(obj), error = function(e) stop("Taxonomy table not found."))
  metadata <- tryCatch(phyloseq::sample_data(obj), error = function(e) stop("Sample metadata not found."))
  otu <- tryCatch(phyloseq::otu_table(obj), error = function(e) stop("OTU/ASV table not found."))

  if (!"spiked.volume" %in% colnames(metadata)) {
    stop("The 'spiked.volume' column is missing in sample metadata.")
  }

  # --- Match spike-ins ---
  matched_taxa <- rownames(tax_data)[tax_data[, "Species"] %in% merged_spiked_species]
  if (length(matched_taxa) == 0 && "Genus" %in% colnames(tax_data)) {
    message("No match found in 'Species' column. Trying 'Genus' column...")
    matched_taxa <- rownames(tax_data)[tax_data[, "Genus"] %in% merged_spiked_species]
  }
  if (length(matched_taxa) == 0) {
    stop("No matching taxa found in 'Species' or 'Genus' columns.")
  }

  spiked_species <- phyloseq::prune_taxa(matched_taxa, obj)
  spiked_species_merged <- phyloseq::merge_taxa(spiked_species, phyloseq::taxa_names(spiked_species))

  # Filter out spike-ins from main object
  filtered_obj <- phyloseq::prune_taxa(!phyloseq::taxa_names(obj) %in% phyloseq::taxa_names(spiked_species), obj)

  # --- Reads ---
  total_reads <- data.frame(
    Sample = phyloseq::sample_names(obj),
    Total_Reads = phyloseq::sample_sums(obj)
  )

  spiked_reads <- data.frame(
    Sample = phyloseq::sample_names(spiked_species_merged),
    Spiked_Reads = phyloseq::sample_sums(spiked_species_merged)
  )

  # --- Scaling factors ---
  scaling_factors <- setNames(rep(1, nrow(total_reads)), total_reads$Sample)
  nonzero <- spiked_reads$Spiked_Reads > 0

  if (any(nonzero)) {
    scaling_factors[nonzero] <- spiked_cells / spiked_reads$Spiked_Reads[nonzero] *
      dplyr::case_when(
        metadata$spiked.volume[nonzero] == 0 ~ 1,
        metadata$spiked.volume[nonzero] == 0.5 ~ 1 / 4,
        metadata$spiked.volume[nonzero] == 1 ~ 1 / 2,
        metadata$spiked.volume[nonzero] == 2 ~ 1,
        metadata$spiked.volume[nonzero] == 3 ~ 1.5,
        metadata$spiked.volume[nonzero] == 4 ~ 2,
        TRUE ~ 1
      )
  }

  # --- Save files ---
  if (!is.null(output_path)) {
    utils::write.csv(total_reads, file.path(output_path, "Total_Reads.csv"), row.names = FALSE)
    utils::write.csv(spiked_reads, file.path(output_path, "Spiked_Reads.csv"), row.names = FALSE)
    utils::write.csv(data.frame(Sample = names(scaling_factors), Scaling_Factor = scaling_factors),
      file.path(output_path, "Scaling_Factors.csv"),
      row.names = FALSE
    )
  }

  # --- Assemble final object ---
  phylo_args <- list(
    phyloseq::otu_table(filtered_obj),
    phyloseq::tax_table(filtered_obj),
    phyloseq::sample_data(filtered_obj)
  )
  if (!is.null(tree)) phylo_args <- c(phylo_args, list(tree))
  if (!is.null(refseq)) phylo_args <- c(phylo_args, list(refseq))

  out_obj <- do.call(phyloseq::phyloseq, phylo_args)

  return(list(
    scaling_factors = scaling_factors,
    filtered_obj = out_obj,
    spiked_species_reads = spiked_reads,
    total_reads = total_reads,
    spiked_species_merged = spiked_species_merged,
    tree = tree,
    refseq = refseq
  ))
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
