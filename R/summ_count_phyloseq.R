#' Summary Statistics of a Phyloseq Object
#'
#' This function computes overall summary statistics (mean, median, standard deviation, standard error, and quantiles) for the OTU table in a phyloseq object.
#'
#' @param physeq A phyloseq object containing the taxonomic and abundance data.
#' @return A data frame with overall summary statistics.
#' @examples
#' # Example usage:
#' # summary_stats <- summ_count_phyloseq(spiked_ITS)
#' # summary_stats <- summ_count_phyloseq(processed_data)
#' @export
summ_count_phyloseq <- function(physeq) {
  # Check if matrixStats package is installed, if not, install it
  if (!require(matrixStats)) {
    install.packages("matrixStats")
  }
  library(matrixStats)
  
  # Extract OTU table
  otu_table <- as.matrix(phyloseq::otu_table(physeq))
  
  # Compute summary statistics
  overall_summary <- data.frame(
    Variable = c("Mean", "Median", "SD", "SE", "Q25", "Q50", "Q75"),
    Value = c(
      mean(rowMeans(otu_table, na.rm = TRUE)),
      median(rowMedians(otu_table, na.rm = TRUE)),
      sd(rowMeans(otu_table, na.rm = TRUE)),
      sd(rowMeans(otu_table, na.rm = TRUE)) / sqrt(nrow(otu_table)),
      quantile(rowMeans(otu_table, na.rm = TRUE), 0.25),
      quantile(rowMeans(otu_table, na.rm = TRUE), 0.50),
      quantile(rowMeans(otu_table, na.rm = TRUE), 0.75)
    )
  )
  
  return(overall_summary)
}


# Example usage:
# summary_stats <- summ_count_phyloseq(spiked_ITS)
# summary_stats <- summ_count_phyloseq(processed_data)
