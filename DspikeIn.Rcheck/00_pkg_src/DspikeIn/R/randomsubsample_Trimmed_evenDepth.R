#' @title Subsampling to an Equal Sequencing Depth
#'
#' @description Performs subsampling to an equal sequencing depth, determined by the sample with
#' the lowest sequencing depth after excluding very low abundant taxa. It rounds down the result
#' to the nearest integer (`floor`). Note: some samples may be lost in this process.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param smalltrim A numeric value specifying the trimming percentage to exclude very low abundant taxa. Default is 0.001.
#' @param replace A logical value indicating whether to sample with replacement. Default is TRUE.
#' @param output_file A character string specifying the output file name for the subsampled object.
#' Default is `NULL`, meaning the object will not be saved.
#'
#' @return A rarefied `phyloseq` or `TreeSummarizedExperiment` object with adjusted sequencing depths.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_ITSOTU", package = "DspikeIn")
#'   tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
#'   rarefied <- randomsubsample_Trimmed_evenDepth(tse_ITSOTU, smalltrim = 0.001)
#'   print(rarefied)
#' }
#'
#' @importFrom phyloseq sample_sums rarefy_even_depth
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @export
randomsubsample_Trimmed_evenDepth <- function(obj,
                                              smalltrim = 0.001,
                                              replace = TRUE,
                                              output_file = NULL) {
  # Validate input
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  }

  # Get OTU table
  otu_matrix <- get_otu_table(obj)

  # Calculate per-sample sequencing depths
  sample_depths <- colSums(otu_matrix)
  sample_depths <- sample_depths[sample_depths > 0]

  if (length(sample_depths) == 0) {
    stop("No valid sample depths found after filtering.")
  }

  # Print sample depth summary
  message("Number of samples: ", length(sample_depths))
  message("Sequencing depth per sample:")
  for (i in seq_along(sample_depths)) {
    sample_name <- names(sample_depths)[i]
    sample_depth <- sample_depths[i]
    message("  ", sample_name, ": ", sample_depth)
  }

  # Sort depths and determine trim index
  sorted_depths <- sort(sample_depths)
  trim_index <- max(1, floor(smalltrim * length(sorted_depths)))
  message("Trim index used: ", trim_index)

  # Compute minimum depth for rarefaction
  samplemin <- sorted_depths[trim_index]
  message("Subsampling to even depth of: ", samplemin)

  # Perform rarefaction
  if (inherits(obj, "phyloseq")) {
    rarefied_obj <- phyloseq::rarefy_even_depth(
      obj,
      sample.size = samplemin,
      rngseed = FALSE,
      replace = replace,
      trimOTUs = TRUE,
      verbose = FALSE
    )
  } else {
    physeq_obj <- convert_tse_to_phyloseq(obj)
    rarefied_physeq <- phyloseq::rarefy_even_depth(
      physeq_obj,
      sample.size = samplemin,
      rngseed = FALSE,
      replace = replace,
      trimOTUs = TRUE,
      verbose = FALSE
    )
    rarefied_obj <- convert_phyloseq_to_tse(rarefied_physeq)
  }

  # Save output if specified
  if (!is.null(output_file)) {
    saveRDS(rarefied_obj, file = output_file)
    message("Rarefied object saved to: ", output_file)
  }

  return(rarefied_obj)
}

# Example usage:
# tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
# ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(tse_ITSOTU, smalltrim = 0.001)
