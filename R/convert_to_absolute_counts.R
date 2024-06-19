#' Convert Relative ASV Counts to Absolute Counts
#'
#' This function converts the relative ASV counts in a phyloseq object to absolute counts by multiplying 
#' the ASV counts by provided scaling factors. The resulting absolute counts are saved as a CSV file.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param scaling_factors A numeric vector of scaling factors to convert relative counts to absolute counts.
#' @param output_dir A character string specifying the directory to save the output file. Default is NULL, which uses the current working directory.
#' @return NULL. The function saves the absolute counts as a CSV file and prints a message indicating the save location.
#' @examples
#' # Example usage:
#' # Assume `spiked_16S` is a phyloseq object and `scaling_factors` is a numeric vector of scaling factors.
#' # convert_to_absolute_counts(spiked_16S, scaling_factors)
#' @export
convert_to_absolute_counts <- function(physeq, scaling_factors, output_dir = NULL) {
  # Load necessary library
  library(phyloseq)
  
  # Change zeros in scaling factors to 1
  scaling_factors[scaling_factors == 0] <- 1
  
  # Convert ASV counts to absolute counts by multiplying ASVs by scaling_factors and rounding the result
  physeq_count <- round(otu_table(physeq) * scaling_factors)
  
  # Check for negative values and replace them with 0
  physeq_count[physeq_count < 0] <- 0
  
  # Set the default output directory to the current working directory if not provided
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Get the file name
  output_file <- file.path(output_dir, "physeq_adj_scaled_AbsoluteCount.csv")
  
  # Save the absolute counts as a CSV file
  write.csv(physeq_count, file = output_file, row.names = TRUE)
  cat("Absolute count data saved to:", output_file, "\n")
}

# Example:
# scaling_factors <- result$scaling_factors
# convert_to_absolute_counts(physeqASV16, scaling_factors)


