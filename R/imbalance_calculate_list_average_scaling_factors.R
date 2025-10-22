#' @title Calculate Per-Sample Scaling Factors for Multiple Spike-in Groups
#'
#' @description
#' Computes per-sample scaling factors for multiple spike-in taxa (e.g.,
#' *Bacillus_spike*, *Flavobacterium_spike*) in either a `phyloseq` or
#' `TreeSummarizedExperiment` object. Handles variable spike-in cell counts per
#' sample and supports `"sum"` or `"max"` OTU merging methods.
#'
#' @details
#' Scaling factors are computed as:
#' \deqn{ScalingFactor = ExpectedSpikeCells / ObservedSpikeReads}
#'
#' For each spike-in group:
#' \enumerate{
#'   \item Identify OTUs matching that spike species via the `Species` column.
#'   \item Merge those OTUs per sample (`sum` or `max`).
#'   \item Divide expected spike cells by observed reads.
#'   \item Average across all spike-in groups to produce one factor per sample.
#' }
#'
#' Uses full matrix preallocation (no incremental vector growth) for Bioconductor
#' compliance. Missing values (zero spike reads) are set to `NA` or `Inf` if
#' `allow_infinite = TRUE`. Samples with all `NA` receive scaling = 1.
#'
#' @param obj A `phyloseq::phyloseq` or
#'   `TreeSummarizedExperiment::TreeSummarizedExperiment` object.
#' @param spiked_species_list A named list of character vectors giving the
#'   spike-in species names (as in `tax_table$Species` or `rowData()`).
#' @param spiked_cells_list A named list (same length as `spiked_species_list`)
#'   containing scalar or named numeric vectors of expected spike-in cells per sample.
#' @param merge_method `"sum"` (default) or `"max"`. Defines how OTUs within a
#'   spike-in group are merged.
#' @param normalize Logical; if `TRUE`, scaling factors are normalized
#'   so that their median equals 1. Default = `TRUE`.
#' @param allow_infinite Logical; if `TRUE`, zero spike reads return `Inf`
#'   instead of `NA`. Default = `FALSE`.
#' @param verbose Logical; if `TRUE`, prints per-group summaries.
#'
#' @return Named numeric vector of scaling factors (one per sample).
#'
#' @importFrom phyloseq otu_table tax_table taxa_are_rows
#' @importFrom SummarizedExperiment assay rowData colnames
#' @importFrom stats median
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   library(phyloseq)
#'
#'   ## Example dataset
#'   otu <- matrix(
#'     c(
#'       6000, 6200, 5900, 6100,
#'       4000, 4200, 3900, 4100,
#'       2000, 1900, 2100, 2050,
#'       1300, 1250, 1350, 1400,
#'        500,  800,  900,  700,   # Flavobacterium_spike
#'        900, 1200, 1100, 1000    # Bacillus_spike
#'     ),
#'     nrow = 6, byrow = TRUE,
#'     dimnames = list(
#'       c("OTU1", "OTU2", "OTU3", "OTU4",
#'         "Flavobacterium_spike", "Bacillus_spike"),
#'       c("S1", "S2", "S3", "S4")
#'     )
#'   )
#'
#'   tax <- data.frame(
#'     Kingdom = rep("Bacteria", 6),
#'     Species = c("OTU1", "OTU2", "OTU3", "OTU4",
#'                 "Flavobacterium_spike", "Bacillus_spike"),
#'     row.names = rownames(otu)
#'   )
#'
#'   #  Fixed: add a column so sample_data is valid
#'   sam <- data.frame(SampleID = c("S1", "S2", "S3", "S4"),
#'                     row.names = c("S1", "S2", "S3", "S4"))
#'
#'   ps <- phyloseq(
#'     otu_table(otu, taxa_are_rows = TRUE),
#'     tax_table(as.matrix(tax)),
#'     sample_data(sam)
#'   )
#'
#'   spiked_species_list <- list(
#'     Flavo = "Flavobacterium_spike",
#'     Bacillus = "Bacillus_spike"
#'   )
#'
#'   spiked_cells_list <- list(
#'     Flavo = c(S1 = 1e7, S2 = 3e7, S3 = 6e7, S4 = 2e7),
#'     Bacillus = c(S1 = 2e7, S2 = 1e7, S3 = 5e7, S4 = 3e7)
#'   )
#'
#'   ## Works for both phyloseq and TSE:
#'   factors_phy <- imbalance_calculate_list_average_scaling_factors(
#'     ps, spiked_species_list, spiked_cells_list, normalize = FALSE
#'   )
#'
#'   tss <- convert_phyloseq_to_tse(ps)
#'   factors_tse <- imbalance_calculate_list_average_scaling_factors(
#'     tss, spiked_species_list, spiked_cells_list, normalize = FALSE
#'   )
#'
#'   all.equal(factors_phy, factors_tse)
#' }
#'
#' @export
imbalance_calculate_list_average_scaling_factors <- function(
    obj,
    spiked_species_list,
    spiked_cells_list,
    merge_method = c("sum", "max"),
    normalize = TRUE,
    allow_infinite = FALSE,
    verbose = FALSE
) {
  merge_method <- match.arg(merge_method)
  
  ## --- Validate input ---
  if (length(spiked_species_list) != length(spiked_cells_list)) {
    stop("'spiked_species_list' and 'spiked_cells_list' must have the same length.")
  }
  
  is_tse <- inherits(obj, "TreeSummarizedExperiment")
  is_phy <- inherits(obj, "phyloseq")
  if (!is_tse && !is_phy)
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  
  ## --- Extract OTU & taxonomy ---
  if (is_tse) {
    otu_mat <- as.matrix(SummarizedExperiment::assay(obj))
    tax_data <- as.data.frame(SummarizedExperiment::rowData(obj))
    sample_names_vec <- SummarizedExperiment::colnames(obj)
  } else {
    otu_mat <- as(phyloseq::otu_table(obj), "matrix")
    ## FIXED: transpose only if taxa_are_rows == FALSE
    if (!phyloseq::taxa_are_rows(obj)) otu_mat <- t(otu_mat)
    tax_data <- as.data.frame(phyloseq::tax_table(obj))
    sample_names_vec <- colnames(otu_mat)
  }
  
  if (!"Species" %in% colnames(tax_data))
    stop("Taxonomy table must contain a 'Species' column.")
  
  n_samples <- length(sample_names_vec)
  n_groups <- length(spiked_species_list)
  
  ## --- Preallocate matrix (efficient) ---
  scaling_factors_matrix <- matrix(
    NA_real_, nrow = n_samples, ncol = n_groups,
    dimnames = list(sample_names_vec, names(spiked_species_list))
  )
  
  ## --- Loop over spike-in groups ---
  for (i in seq_len(n_groups)) {
    spike_group <- names(spiked_species_list)[i]
    spiked_species <- spiked_species_list[[i]]
    expected_cells <- spiked_cells_list[[i]]
    
    matched_otus <- which(tax_data$Species %in% spiked_species)
    if (!length(matched_otus)) {
      warning("No OTUs matched for spike group: ", spike_group)
      next
    }
    
    ## Different indexing for TSE vs phyloseq
    spike_abund <- if (is_tse) {
      otu_mat[matched_otus, , drop = FALSE]
    } else {
      otu_mat[matched_otus, , drop = FALSE]
    }
    
    merged_abundance <- switch(
      merge_method,
      sum = colSums(spike_abund, na.rm = TRUE),
      max = apply(spike_abund, 2, max, na.rm = TRUE)
    )
    
    expected_vector <- if (length(expected_cells) == 1) {
      rep(expected_cells, n_samples)
    } else if (all(names(expected_cells) %in% sample_names_vec)) {
      expected_cells[sample_names_vec]
    } else {
      stop("Each expected_cells must be scalar or named vector matching sample names.")
    }
    
    scaling_vec <- expected_vector / merged_abundance
    zero_or_na <- merged_abundance == 0 | is.na(merged_abundance)
    scaling_vec[zero_or_na] <- if (allow_infinite) Inf else NA_real_
    
    scaling_factors_matrix[, i] <- scaling_vec
    
    if (verbose) {
      cat("\n Spike group:", spike_group,
          "\nMatched taxa:", paste(spiked_species, collapse = ", "),
          "\nMean reads:", round(mean(merged_abundance, na.rm = TRUE), 2),
          "\nMean scaling factor:", round(mean(scaling_vec, na.rm = TRUE), 2), "\n")
    }
  }
  
  ## --- Aggregate across groups ---
  avg_scaling <- rowMeans(scaling_factors_matrix, na.rm = TRUE)
  avg_scaling[is.na(avg_scaling)] <- 1
  
  ## --- Normalize safely ---
  if (normalize) {
    finite_vals <- avg_scaling[is.finite(avg_scaling)]
    med_val <- if (length(finite_vals)) stats::median(finite_vals) else 1
    avg_scaling <- avg_scaling / med_val
  }
  
  ## --- Verbose summary ---
  if (verbose) {
    inf_count <- sum(is.infinite(avg_scaling))
    if (inf_count > 0) {
      cat(inf_count, "sample(s) had zero spike reads-scaling =",
          ifelse(allow_infinite, "Inf", "NA"), "\n")
    }
    cat("Calculated scaling factors for", n_samples,
        "samples using", n_groups, "spike-in groups (merge =", merge_method, ").\n")
  }
  
  return(round(avg_scaling, 6))
}




