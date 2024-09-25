#' Subsampling to an Equal Sequencing Depth
#'
#' This function performs subsampling to an equal sequencing depth, determined by the sample with the lowest sequencing depth after excluding very low abundant taxa.
#' It rounds down the result to the nearest integer (floor). Caution: you might lose samples.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param smalltrim A numeric value specifying the trimming percentage to exclude very low abundant taxa. Default is 0.001.
#' @param replace A logical value indicating whether to sample with replacement. Default is TRUE.
#' @param output_file A character string specifying the output file name for the subsampled phyloseq object. Default is "randomsubsample_Trimmed_evenDepth.rds".
#' @return A phyloseq object with subsampled sequencing depths.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Perform subsampling
#'   spiked_ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(spiked_ITS, smalltrim = 0.001)
#' }
#' }
#' @importFrom phyloseq sample_sums rarefy_even_depth
#' @export
randomsubsample_Trimmed_evenDepth <- function(physeq, smalltrim = 0.001, replace = TRUE, output_file = "randomsubsample_Trimmed_evenDepth.rds") {
  suppressMessages({
    # Calculate the sample sizes and remove zeros
    sample_depths <- sort(phyloseq::sample_sums(physeq))
    sample_depths <- sample_depths[sample_depths > 0]
    
    # Ensure we have valid sample sizes
    if (length(sample_depths) == 0) {
      stop("No valid sample depths found after removing zero counts. Check the phyloseq object.")
    }
    
    # Calculate the trim index
    trim_index <- max(1, floor(smalltrim * length(sample_depths)))
    
    # Debugging 
    cat("Number of samples:", length(sample_depths), "\n")
    cat("Sample depths (sorted):", sample_depths, "\n")
    cat("Trim index:", trim_index, "\n")
    
    # If out of bounds
    if (trim_index >= length(sample_depths)) {
      stop("Calculated trim_index is out of bounds. Adjust the smalltrim value.")
    }
    
    # Calculate samplemin
    samplemin <- sample_depths[trim_index]
    
    # Debugging information
    cat("Calculated samplemin:", samplemin, "\n")
    
    # Ensure samplemin is greater than zero
    if (samplemin <= 0) {
      stop("Calculated samplemin is less than or equal to zero. Adjust the smalltrim value.")
    }
    
    # Rarefy to even depth
    physeq_evenDepth <- phyloseq::rarefy_even_depth(physeq, samplemin, rngseed = FALSE, replace = replace, trimOTUs = TRUE)
    
    saveRDS(physeq_evenDepth, file = output_file)
    cat("Rarefied phyloseq object saved to:", output_file, "\n")
    
    return(physeq_evenDepth)
  })
}

# Example usage:
# spiked_ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(spiked_ITS, smalltrim = 0.001)
