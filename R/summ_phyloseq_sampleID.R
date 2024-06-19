#' Generate Summary Statistics for Each Sample in a Phyloseq Object
#'
#' This function calculates summary statistics (mean, median, standard deviation, standard error, and quartiles) for each sample in a phyloseq object.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @return A data frame containing the summary statistics for each sample, with columns for Sample_ID, mean, median, standard deviation, standard error, and quartiles (Q25, Q50, Q75).
#' @examples
#' # Example usage:
#' # Generate summary statistics for a phyloseq object
#' # summary_stats <- summ_phyloseq_sampleID(spiked_16S)
#' @export
summ_phyloseq_sampleID <- function(physeq) {
  otu_table <- otu_table(physeq)
  summary_stats <- data.frame(
    Sample_ID = colnames(otu_table),
    Mean = apply(otu_table, 2, mean),
    Median = apply(otu_table, 2, median),
    SD = apply(otu_table, 2, sd),
    SE = apply(otu_table, 2, function(x) sd(x) / sqrt(length(x))),
    Q25 = apply(otu_table, 2, quantile, probs = 0.25),
    Q50 = apply(otu_table, 2, quantile, probs = 0.5),
    Q75 = apply(otu_table, 2, quantile, probs = 0.75)
  )
  return(summary_stats)
}

# Example usage:
# Generate summary statistics for a phyloseq object
#summary_stats <- summ_phyloseq_sampleID(spiked_16S)
