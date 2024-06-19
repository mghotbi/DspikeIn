#' Summarize Phyloseq Data Based on ASV_ID
#'
#' This function generates summary statistics (mean, median, standard deviation, standard error, and quantiles)
#' for each ASV (Amplicon Sequence Variant) in a phyloseq object.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @return A data frame containing summary statistics for each ASV.
#' @examples
#' # Summarize the phyloseq data based on ASV_ID
#' summary_stats <- summ_phyloseq_ASV_OTUID(spiked_16S)
#' print(summary_stats)
#' @export
summ_phyloseq_ASV_OTUID <- function(physeq) {
  # Extract the OTU table from the phyloseq object
  otu_table <- otu_table(physeq)
  
  # Calculate summary statistics for each ASV
  summary_stats <- data.frame(
    ASV_ID = rownames(otu_table),          # ASV identifier
    Mean = apply(otu_table, 1, mean),      # Mean abundance
    Median = apply(otu_table, 1, median),  # Median abundance
    SD = apply(otu_table, 1, sd),          # Standard deviation
    SE = apply(otu_table, 1, function(x) sd(x) / sqrt(length(x))),  # Standard error
    Q25 = apply(otu_table, 1, quantile, probs = 0.25),  # 25th percentile (1st quartile)
    Q50 = apply(otu_table, 1, quantile, probs = 0.5),   # 50th percentile (2nd quartile / median)
    Q75 = apply(otu_table, 1, quantile, probs = 0.75)   # 75th percentile (3rd quartile)
  )
  
  # Return the summary statistics data frame
  return(summary_stats)
}

# Example usage:
# Summarize the phyloseq data based on ASV_ID
#summary_stats <- summ_phyloseq_ASV_OTUID(spiked_16S)
#print(summary_stats)
