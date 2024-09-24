#' Summary Statistics of a Phyloseq Object
#'
#' This function computes overall summary statistics (mean, median, standard deviation, standard error, and quantiles) for the OTU table in a phyloseq object.
#'
#' @param physeq A phyloseq object containing the taxonomic and abundance data.
#' @return A data frame with overall summary statistics.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   summary_stats <- summ_count_phyloseq(spiked_ITS)
#'   summary_stats <- summ_count_phyloseq(processed_data)
#' }
#' }
#' @importFrom matrixStats rowMeans2 rowMedians
#' @importFrom stats median sd quantile
#' @export
summ_count_phyloseq <- function(physeq) {
  suppressMessages({
    # Extract OTU table
    otu_table <- as.matrix(phyloseq::otu_table(physeq))
    
    # Compute summary statistics
    overall_summary <- data.frame(
      Variable = c("Mean", "Median", "SD", "SE", "Q25", "Q50", "Q75"),
      Value = c(
        base::mean(matrixStats::rowMeans2(otu_table, na.rm = TRUE)),
        stats::median(matrixStats::rowMedians(otu_table, na.rm = TRUE)),
        stats::sd(matrixStats::rowMeans2(otu_table, na.rm = TRUE)),
        stats::sd(matrixStats::rowMeans2(otu_table, na.rm = TRUE)) / base::sqrt(nrow(otu_table)),
        stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.25),
        stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.50),
        stats::quantile(matrixStats::rowMeans2(otu_table, na.rm = TRUE), 0.75)
      )
    )
    
    return(overall_summary)
  })
}

# Example usage:
# summary_stats <- summ_count_phyloseq(spiked_ITS)
# summary_stats <- summ_count_phyloseq(processed_data)
