#' @title Calculate Average Scaling Factors for Multiple Spiked Species
#'
#' @description This function calculates scaling factors for multiple spiked species in a \code{phyloseq} or
#' \code{TreeSummarizedExperiment} (TSE) object. It merges ASVs/OTUs for each species if necessary,
#' averages the scaling factors, and returns the averaged scaling factors for each OTU.
#' Different spiked cell counts can be provided for each set of spiked species.
#' If an OTU is not associated with any spiked species, a default scaling factor (1) is assigned.
#' Scaling factors are rounded to the specified number of decimal places.
#'
#' @param obj A \code{phyloseq} or \code{TreeSummarizedExperiment} (TSE) object containing microbial abundance data.
#' @param spiked_species_list A list of character vectors. Each vector contains the spiked species
#' (by taxon names) to be merged for calculating scaling factors for each group.
#' @param spiked_cells_list A numeric vector specifying the number of spiked cells corresponding to each group in `spiked_species_list`.
#' @param merge_method A character string specifying how to merge ASVs/OTUs for each group of spiked species.
#' Accepted values are \code{"sum"} or \code{"max"}. Default is \code{"sum"}.
#' @return A numeric vector of averaged and rounded scaling factors for each OTU in the object.
#' OTUs not associated with any spiked species are assigned a default scaling factor of 1.
#'
#' @importFrom phyloseq taxa_names otu_table tax_table
#' @importFrom SummarizedExperiment assay rowData
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq", package = "DspikeIn")
#'
#'   # Step 1: Define the spiked species list and corresponding cell counts
#'   spiked_species_list <- list(
#'     c("Pseudomonas aeruginosa"),
#'     c("Escherichia coli"),
#'     c("Clostridium difficile")
#'   )
#'
#'   spiked_cells_list <- c(10000, 20000, 15000)
#'
#'   # Step 2: Apply the function to a phyloseq object
#'   scaling_factors <- calculate_list_average_scaling_factors(
#'     physeq,
#'     spiked_species_list,
#'     spiked_cells_list,
#'     merge_method = "sum"
#'   )
#'
#'   # Step 3: Print the results
#'   print(scaling_factors)
#' }
#' }
#' @export
calculate_list_average_scaling_factors <- function(obj, spiked_species_list, spiked_cells_list, merge_method = c("sum", "max")) {

  # Ensure correct merge_method input
  merge_method <- match.arg(merge_method)

  # Validate input lengths
  if (length(spiked_species_list) != length(spiked_cells_list)) {
    stop("\U0000274C The length of spiked_species_list must match the length of spiked_cells_list.")
  }

  # Determine object type and use appropriate accessor functions
  otu_table <- if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::assay(obj)
  } else {
    phyloseq::otu_table(obj)
  }

  tax_data <- if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::rowData(obj)
  } else {
    phyloseq::tax_table(obj)
  }

  otu_names <- rownames(otu_table)

  # Initialize scaling factor storage
  scaling_factors_total <- rep(0, length(otu_names))
  count_contributions <- rep(0, length(otu_names))

  # Loop through each group of spiked species and calculate scaling factors
  for (i in seq_along(spiked_species_list)) {

    # Identify OTUs corresponding to the spiked species
    spiked_species <- spiked_species_list[[i]]
    matched_otus <- which(tax_data[, "Species"] %in% spiked_species)

    # Ensure OTUs were found
    if (length(matched_otus) == 0) {
      warning(paste("No OTUs matched for spiked species:", paste(spiked_species, collapse = ", ")))
      next
    }

    # Extract relevant abundance data
    spiked_abundances <- otu_table[matched_otus, , drop = FALSE]

    # Merge ASVs/OTUs based on specified method
    merged_abundance <- if (merge_method == "sum") {
      rowSums(spiked_abundances, na.rm = TRUE)
    } else {
      apply(spiked_abundances, 2, max, na.rm = TRUE)
    }

    # Calculate total observed abundance
    total_abundance_spiked <- sum(merged_abundance)

    # Avoid division by zero
    if (total_abundance_spiked == 0) {
      warning(paste("Zero total abundance for spiked species:", paste(spiked_species, collapse = ", ")))
      next
    }

    # Compute scaling factor
    scaling_factor <- spiked_cells_list[i] / total_abundance_spiked

    # Assign scaling factors to matched OTUs
    scaling_factors_total[matched_otus] <- scaling_factors_total[matched_otus] + scaling_factor
    count_contributions[matched_otus] <- count_contributions[matched_otus] + 1
  }

  # Compute average scaling factors
  average_scaling_factors <- scaling_factors_total / count_contributions

  # Replace NAs with default scaling factor (1)
  average_scaling_factors[is.na(average_scaling_factors)] <- 1

  # Round values for consistency
  average_scaling_factors <- round(average_scaling_factors, digits = 2)

  return(average_scaling_factors)
}

# Example usage:
# Step 1: Define the spiked species list and corresponding cell counts
#spiked_species_list <- list(
#  c("Pseudomonas aeruginosa"),
#  c("Escherichia coli"),
#  c("Clostridium difficile")
#)

# spiked_cells_list <- c(10000, 20000, 15000)

# Step 2: Apply the function to a phyloseq object
#scaling_factors <- calculate_list_average_scaling_factors(
# physeq,
# spiked_species_list,
#  spiked_cells_list,
#  merge_method = "sum")

# Step 3: Print the results
#print(scaling_factors)

#
# # Step 4: build the phyloseq
# otu_table_ps <- phyloseq::otu_table(otu_data, taxa_are_rows = TRUE)
# tax_table_ps <- phyloseq::tax_table(taxa_matrix)
# sample_data_ps <- phyloseq::sample_data(sample_data)
#
# physeq <- phyloseq(otu_table_ps, tax_table_ps, sample_data_ps)
#
# # tidy up
# physeq<- tidy_phyloseq(physeq)
#
# spiked_species_list <- list(
#   c("Pseudomonas aeruginosa"),
#   c("Escherichia coli"),
#   c("Clostridium difficile")
# )
#
# spiked_cells_list <- c(10000, 20000, 15000)
#
# # Step 6: Calculate the scaling factors after merging the redundant spikein species
# scaling_factors <- calculate_list_average_scaling_factors(merged_physeq_sum,
# spiked_species_list, spiked_cells_list, merge_method = "sum") # or max
# # Print the scaling factors for each OTU
# print(scaling_factors)
