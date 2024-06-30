#' Filter Taxa from a Phyloseq Object Based on Custom Thresholds
#'
#' This function filters taxa from a phyloseq object based on custom thresholds for percentage of samples, mean abundance, count, and relative abundance.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param threshold_percentage A numeric value specifying the minimum percentage of samples in which a taxon must be present to be retained. Default is 0.5.
#' @param threshold_mean_abundance A numeric value specifying the minimum mean abundance of a taxon to be retained. Default is 0.001.
#' @param threshold_count A numeric value specifying the minimum count of a taxon in a sample to be considered present. Default is 10.
#' @param threshold_relative_abundance A numeric value specifying the minimum relative abundance of a taxon to be retained. Default is NULL.
#' @return A phyloseq object containing only the taxa that meet the specified thresholds.
#' @examples
#' # Example usage with custom thresholds
#' # FT <- relativized_filtered_taxa(spiked_16S, threshold_percentage = 0.6, 
#' # threshold_mean_abundance = 0.0005,
#' # threshold_count = 5, threshold_relative_abundance = 0.01)
#' @export
relativized_filtered_taxa <- function(physeq, 
                                      threshold_percentage = 0.5, 
                                      threshold_mean_abundance = 0.001, 
                                      threshold_count = 10,
                                      threshold_relative_abundance = NULL) {
  
  nsamples <- nsamples(physeq)
  sample_sum <- sample_sums(physeq)
  
  filter_function <- function(x) {
    (sum(x > threshold_count) > nsamples * threshold_percentage) | 
      ((sum(x > threshold_count) > (nsamples * 0.1)) & (mean(x / sample_sum) > threshold_mean_abundance) & (max(x / sample_sum) > threshold_relative_abundance))
  }
  
  two_way_filtered <- filter_taxa(physeq, filter_function, prune = TRUE)
  return(two_way_filtered)
}

# Example usage with custom thresholds
# FT <- relativized_filtered_taxa(spiked_16S, threshold_percentage = 0.6, 
#threshold_mean_abundance = 0.0005,threshold_count = 5, threshold_relative_abundance = 0.01)
