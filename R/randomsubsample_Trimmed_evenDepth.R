#subsampling to an equal sequencing depth, determined by the sample with the lowest sequencing depth after excluding very low abundant taxa
#rounds down the result to the nearest integer/floor. Caution you might loose samples
randomsubsample_Trimmed_evenDepth <- function(physeq, smalltrim = 0.001, replace = TRUE, output_file = "randomsubsample_Trimmed_evenDepth.rds") {
  # Calculate the sample sizes and remove zeros
  sample_depths <- sort(sample_sums(physeq))
  sample_depths <- sample_depths[sample_depths > 0]
  
  # Ensure we have valid sample sizes
  if (length(sample_depths) == 0) {
    stop("No valid sample depths found after removing zero counts. Check the phyloseq obj.")
  }
  
  # Calculate the trim index
  trim_index <- max(1, floor(smalltrim * length(sample_depths)))
  
  # Debugging 
  cat("Number of samples:", length(sample_depths), "\n")
  cat("Sample depths (sorted):", sample_depths, "\n")
  cat("Trim index:", trim_index, "\n")
  
  # If out of bounds?
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
  physeq_evenDepth <- rarefy_even_depth(physeq, samplemin, rngseed = FALSE, replace = replace, trimOTUs = TRUE)
  
  saveRDS(physeq_evenDepth, file = output_file)
  cat("Rarefied phyloseq object saved to:", output_file, "\n")
  
  return(physeq_evenDepth)
}

# Example usage
#spiked_ITS_evenDepth <- randomsubsample_Trimmed_evenDepth(spiked_ITS, smalltrim = 0.001)


