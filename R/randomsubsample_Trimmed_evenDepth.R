#' @title Subsampling to an Equal Sequencing Depth
#'
#' @description Performs subsampling to an equal sequencing depth, determined by the sample with
#' the lowest sequencing depth after excluding very low abundant taxa. It rounds down the result
#' to the nearest integer (floor). Caution: you might lose samples.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param smalltrim A numeric value specifying the trimming percentage to exclude very low abundant taxa. Default is 0.001.
#' @param replace A logical value indicating whether to sample with replacement. Default is TRUE.
#' @param output_file A character string specifying the output file name for the subsampled object. Default is "randomsubsample_Trimmed_evenDepth.rds".
#' @return A `phyloseq` or `TreeSummarizedExperiment` object with subsampled sequencing depths.
#'
#' @importFrom phyloseq sample_sums rarefy_even_depth
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_ITSOTU", package = "DspikeIn")
#'
#'   # Convert phyloseq to TreeSummarizedExperiment (TSE)
#'   tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
#'
#'   # Perform random subsampling with trimmed even depth
#'   ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(tse_ITSOTU, smalltrim = 0.001)
#' }
#' }
#' @export
randomsubsample_Trimmed_evenDepth <- function(obj, smalltrim = 0.001, replace = TRUE, output_file = "randomsubsample_Trimmed_evenDepth.rds") {
  suppressMessages({

    #  Detect input type and extract OTU table
    otu_matrix <- get_otu_table(obj)

    #  Detect sample data format (to maintain metadata integrity)
    sample_data_df <- get_sample_data(obj)

    #  Calculate the sample sizes and remove zeros
    sample_depths <- sort(colSums(otu_matrix))  # Use column sums for sample depths
    sample_depths <- sample_depths[sample_depths > 0]

    #  Ensure valid sample sizes exist
    if (length(sample_depths) == 0) {
      stop("	\U0000274C No valid sample depths found after removing zero counts. Check the input object.")
    }

    #  Calculate the trim index
    trim_index <- max(1, floor(smalltrim * length(sample_depths)))

    # 🛠 Debugging information
    cat("Number of samples:", length(sample_depths), "\n")
    cat("Sample depths (sorted):", sample_depths, "\n")
    cat("Trim index:", trim_index, "\n")

    #  Ensure trim index is within bounds
    if (trim_index >= length(sample_depths)) {
      stop("	\U0000274C Calculated trim_index is out of bounds. Adjust the smalltrim value.")
    }

    #  Compute minimum depth for rarefaction
    samplemin <- sample_depths[trim_index]

    #  Debugging information
    cat("	\U0001F9EE Calculated samplemin:", samplemin, "\n")

    #  Ensure `samplemin` is valid
    if (samplemin <= 0) {
      stop("\U0000274C Calculated samplemin is less than or equal to zero. Adjust the smalltrim value.")
    }

    # ️ Rarefy to even depth using phyloseq function
    if (inherits(obj, "phyloseq")) {
      obj_rarefied <- phyloseq::rarefy_even_depth(obj, samplemin, rngseed = FALSE, replace = replace, trimOTUs = TRUE)
    } else if (inherits(obj, "TreeSummarizedExperiment")) {
      # Convert TSE → phyloseq, rarefy, then convert back
      obj_physeq <- convert_tse_to_phyloseq(obj)
      obj_rarefied <- phyloseq::rarefy_even_depth(obj_physeq, samplemin, rngseed = FALSE, replace = replace, trimOTUs = TRUE)
      obj_rarefied <- convert_phyloseq_to_tse(obj_rarefied)
    } else {
      stop("	\U0000274C Unsupported object type: must be 'phyloseq' or 'TreeSummarizedExperiment'.")
    }

    #  Save output
    saveRDS(obj_rarefied, file = output_file)
    cat("\U0001F4C2 Rarefied object saved to:", output_file, "\n")

    return(obj_rarefied)
  })
}


# Example usage:
# tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
# ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(tse_ITSOTU, smalltrim = 0.001)
