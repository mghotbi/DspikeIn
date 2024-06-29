#' Convert Relative ASV Counts to Absolute Counts
#'
#' This function converts the relative ASV counts in a phyloseq object to absolute counts by multiplying 
#' the ASV counts by provided scaling factors. The resulting absolute counts are saved as a CSV file and 
#' also returned as a list containing both the data frame of absolute counts and the modified phyloseq object. 
#' Taxa with NA or empty values in the taxonomic table are removed.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param scaling_factors A numeric vector of scaling factors to convert relative counts to absolute counts.
#' @param output_dir A character string specifying the directory to save the output file. Default is NULL, which uses the current working directory.
#' @return A list containing the data frame of absolute counts and the modified phyloseq object.
#' @examples
#' # Example usage:
#' # Assume `spiked_16S` is a phyloseq object and `scaling_factors` is a numeric vector of scaling factors.
#' # absolute <- convert_to_absolute_counts(spiked_16S, scaling_factors)
#' # absolute_counts <- absolute$absolute_counts
#' # physeq_obj <- absolute$physeq_obj
#' @export
convert_to_absolute_counts <- function(physeq, scaling_factors, output_dir = NULL) {
  # Load necessary library
  library(phyloseq)
  
  # Remove taxa with NA or empty values in the taxonomic table
  physeq <- subset_taxa(physeq, apply(tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
  
  # Change zeros in scaling factors to 1
  scaling_factors[scaling_factors == 0] <- 1
  
  # Convert ASV counts to absolute counts by multiplying ASVs by scaling_factors and rounding the result
  physeq_count <- round(otu_table(physeq) * scaling_factors)
  
  # Check for negative or NA values and replace them with 0
  physeq_count[physeq_count < 0 | is.na(physeq_count)] <- 0
  
  # Set the OTU table in the phyloseq object with the absolute counts
  otu_table(physeq) <- otu_table(physeq_count)
  
  # Set the default output directory to the current working directory if not provided
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Get the file name
  output_file <- file.path(output_dir, "physeq_adj_scaled_AbsoluteCount.csv")
  
  # Save the absolute counts as a CSV file
  write.csv(otu_table(physeq), file = output_file, row.names = TRUE)
  cat("Absolute count data saved to:", output_file, "\n")
  
  # Return a list containing the data frame of absolute counts and the modified phyloseq object
  return(list(
    absolute_counts = as.data.frame(otu_table(physeq)),
    physeq_obj = physeq
  ))
}

# Example:
# scaling_factors <- absolute$scaling_factors
#absolute <- convert_to_absolute_counts(physeq_16SASV, scaling_factors)
#absolute_counts <- absolute$absolute_counts
#physeq_absolute_abundance <- absolute$physeq_obj
