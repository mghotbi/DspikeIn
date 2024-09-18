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
#' @importFrom phyloseq subset_taxa tax_table otu_table
#' @importFrom utils write.csv
#' @examples
#' \dontrun{
#' # Example usage:
#' #' # Assume `merged_physeq_sum` is a phyloseq object with one ASV/OTU 
#' # resulting from spiked species, and `scaling_factors` is a 
#' # numeric vector of scaling factors.
#' merged_spiked_species <- c("Tetragenococcus_halophilus")
#' 
#' # Calculate scaling factors and generate a report
#' result <- calculate_spikeIn_factors(merged_physeq_sum, 1874, merged_spiked_species)
#' 
#' # Convert the phyloseq object to absolute counts using scaling factors
#' absolute <- convert_to_absolute_counts(merged_physeq_sum, scaling_factors)
#' absolute_counts <- absolute$absolute_counts
#' physeq_obj <- absolute$physeq_obj
#' 
#' # The 'absolute_counts' variable contains the converted counts, and 'physeq_obj' 
#' #is the updated phyloseq object
#' }
#' @export
convert_to_absolute_counts <- function(physeq, scaling_factors, output_dir = NULL) {
  suppressMessages({
    # Ensure necessary libraries are available
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("Package 'phyloseq' is required but not installed.")
    }
    
    # Remove taxa with NA or empty values in the taxonomic table
    physeq <- phyloseq::subset_taxa(physeq, apply(phyloseq::tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
    
    # Change zeros in scaling factors to 1
    scaling_factors[scaling_factors == 0] <- 1
    
    # Convert ASV counts to absolute counts by multiplying ASVs by scaling_factors and rounding the result
    physeq_count <- round(phyloseq::otu_table(physeq) * scaling_factors)
    
    # Check for negative or NA values and replace them with 0
    physeq_count[physeq_count < 0 | is.na(physeq_count)] <- 0
    
    # Set the OTU table in the phyloseq object with the absolute counts
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(physeq_count)
    
    # Set the default output directory to the current working directory if not provided
    if (is.null(output_dir)) {
      output_dir <- getwd()
    }
    
    # Get the file name
    output_file <- file.path(output_dir, "physeq_adj_scaled_AbsoluteCount.csv")
    
    # Save the absolute counts as a CSV file
    utils::write.csv(phyloseq::otu_table(physeq), file = output_file, row.names = TRUE)
    cat("Absolute count data saved to:", output_file, "\n")
    
    # Return a list containing the data frame of absolute counts and the modified phyloseq object
    return(list(
      absolute_counts = as.data.frame(phyloseq::otu_table(physeq)),
      physeq_obj = physeq
    ))
  })
}

# Example:
# merged_spiked_species <- c("Tetragenococcus_halophilus")
# # Calculate scaling factors and generate the report
# result <- calculate_spikeIn_factors(merged_physeq_sum, 1874, merged_spiked_species)
# scaling_factors <- result$scaling_factors
# absolute <- convert_to_absolute_counts(merged_physeq_sum, scaling_factors)
# absolute_counts <- absolute$absolute_counts
# physeq_absolute <- absolute$physeq_obj
