#' @title Calculate Sample-specific Average Scaling Factors for Multiple Spike-in Groups
#'
#' @description
#' Computes sample-specific scaling factors for multiple groups of spiked species
#' in a `phyloseq` or `TreeSummarizedExperiment` object. Each group can have its own
#' expected spike-in cell count. Scaling factors are calculated per sample and averaged across groups.
#' Missing spike-in observations in a sample will be handled gracefully by averaging available groups.
#'
#' @details
#' The function assumes that the taxonomy table has a `Species` column.
#' The output is suitable for downstream absolute quantification pipelines.
#' OTUs belonging to each spike-in group will be merged using the specified `merge_method`
#' ("sum" or "max") to obtain a group-specific spike-in abundance in each sample.
#'
#' If a sample does not contain any spike-in sequences, a scaling factor of 1 is assigned.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param spiked_species_list A list of character vectors. Each vector contains taxon names (at species level) for one spike-in group.
#' @param spiked_cells_list A numeric vector specifying the expected number of spike-in cells for each group.
#' The order must match `spiked_species_list`.
#' @param merge_method Character. Either `"sum"` or `"max"`. Controls how OTUs of each spike-in group are merged.
#' @return A named numeric vector of sample-specific scaling factors.
#'
#' @section Notes:
#' - This function does not modify the input object.
#' - The returned scaling factors are intended to be used for absolute abundance normalization.
#'
#' @importFrom phyloseq taxa_names otu_table tax_table
#' @importFrom SummarizedExperiment assay rowData
#' @export
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq", package = "DspikeIn")
#'
#'   # Step 1: Define the spike-in species groups and associated cell counts
#'   spiked_species_list <- list(
#'     c("Pseudomonas aeruginosa"),
#'     c("Escherichia coli"),
#'     c("Clostridium difficile")
#'   )
#'
#'   spiked_cells_list <- c(10000, 20000, 15000)
#'
#'   # Step 2: Compute the scaling factors
#'   scaling_factors <- calculate_list_average_scaling_factors(
#'     physeq,
#'     spiked_species_list,
#'     spiked_cells_list,
#'     merge_method = "sum"
#'   )
#'
#'   # Step 3: Inspect scaling factors
#'   print(scaling_factors)
#' }
#' @export
calculate_list_average_scaling_factors <- function(obj, spiked_species_list, spiked_cells_list, merge_method = c("sum", "max")) {
  merge_method <- match.arg(merge_method)

  # Validate inputs
  if (length(spiked_species_list) != length(spiked_cells_list)) {
    stop("'spiked_species_list' and 'spiked_cells_list' must have the same length.")
  }

  # Extract OTU and taxonomy tables
  is_tse <- inherits(obj, "TreeSummarizedExperiment")

  otu_mat <- if (is_tse) {
    SummarizedExperiment::assay(obj)
  } else {
    as(phyloseq::otu_table(obj), "matrix")
  }

  tax_data <- if (is_tse) {
    as.data.frame(SummarizedExperiment::rowData(obj))
  } else {
    as.data.frame(phyloseq::tax_table(obj))
  }

  sample_names_vec <- colnames(otu_mat)
  n_samples <- ncol(otu_mat)

  # Initialize matrices
  scaling_factors_matrix <- matrix(NA, nrow = n_samples, ncol = length(spiked_species_list))
  rownames(scaling_factors_matrix) <- sample_names_vec

  # Loop over each spike-in group
  for (i in seq_along(spiked_species_list)) {
    spiked_species <- spiked_species_list[[i]]
    expected_cells <- spiked_cells_list[i]

    # Identify matched OTUs by Species
    matched_otus <- which(tax_data$Species %in% spiked_species)

    if (length(matched_otus) == 0) {
      warning(sprintf("No OTUs matched for: %s", paste(spiked_species, collapse = ", ")))
      next
    }

    spikein_abund <- otu_mat[matched_otus, , drop = FALSE]

    # Merge OTUs per sample
    merged_abundance <- if (merge_method == "sum") {
      colSums(spikein_abund, na.rm = TRUE)
    } else {
      apply(spikein_abund, 2, max, na.rm = TRUE)
    }

    # Avoid divide by zero
    scaling_vec <- ifelse(merged_abundance > 0,
      expected_cells / merged_abundance,
      NA
    )

    scaling_factors_matrix[, i] <- scaling_vec
  }

  # Calculate average scaling factor per sample (excluding NAs)
  average_scaling <- rowMeans(scaling_factors_matrix, na.rm = TRUE)

  # Replace any NA (i.e., sample had no spike-in from any group) with default 1
  average_scaling[is.na(average_scaling)] <- 1

  # Round for consistency
  average_scaling <- round(average_scaling, digits = 4)

  names(average_scaling) <- sample_names_vec
  return(average_scaling)
}



# Example usage:
# Step 1: Define the spiked species list and corresponding cell counts
# spiked_species_list <- list(
#  c("Pseudomonas aeruginosa"),
#  c("Escherichia coli"),
#  c("Clostridium difficile")
# )

# spiked_cells_list <- c(10000, 20000, 15000)

# Step 2: Apply the function to a phyloseq object
# scaling_factors <- calculate_list_average_scaling_factors(
# physeq,
# spiked_species_list,
#  spiked_cells_list,
#  merge_method = "sum")

# Step 3: Print the results
# print(scaling_factors)

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
# 
