#' Adjust Prevalence in a Phyloseq Object
#'
#' This function adjusts the prevalence of ASVs (Amplicon Sequence Variants) in a phyloseq object
#' based on a specified method to determine the number of reads (min, mean, median, or max).
#' It then prunes the ASVs that do not meet the prevalence threshold and saves the adjusted phyloseq object.
#'
#' @param spiked_16S A phyloseq object containing the taxonomic and abundance data.
#' @param method A character string specifying the method to calculate the number of reads. Options are "min", "mean", "median", and "max". Default is "min".
#' @param output_file A character string specifying the file path to save the adjusted phyloseq object. Default is "adjusted_prevalence_physeq.rds".
#' @return A phyloseq object with adjusted prevalence.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Adjust prevalence using the "min" method
#'   spiked_16S_min <- adjusted_prevalence(physeq_16S_adj_scaled, method = "min")
#'   
#'   # Adjust prevalence using the "max" method
#'   spiked_16S_max <- adjusted_prevalence(physeq_16S_adj_scaled, method = "max")
#' }
#' }
#' @export
adjusted_prevalence <- function(spiked_16S, method = "min", output_file = "adjusted_prevalence_physeq.rds") {
  suppressMessages({
    # Define the method to calculate the number of reads
    method <- tolower(method)
    
    # Find the number of reads based on the specified method
    reads <- switch(method,
                    "min" = min(phyloseq::taxa_sums(spiked_16S)),
                    "mean" = mean(phyloseq::taxa_sums(spiked_16S)),
                    "median" = median(phyloseq::taxa_sums(spiked_16S)),
                    "max" = max(phyloseq::taxa_sums(spiked_16S)),
                    stop("Invalid method. Please choose 'min', 'mean', 'median', or 'max'."))
    
    # Number of reads
    cat("Number of reads chosen by method", method, ":", reads, "\n")
    
    # Calculate the total number of reads for each sample
    sampsums <- phyloseq::sample_sums(spiked_16S)
    
    # Trim ASVs that do not appear in many samples (prevalence)
    samobs <- apply(phyloseq::otu_table(spiked_16S), 1, function(x) sum(x > 0))  # Count ASV occurrences
    
    # Create a dataframe with prevalence and sums
    otudf <- data.frame(prev = samobs, sums = phyloseq::taxa_sums(spiked_16S))
    otudf <- otudf[order(-otudf$prev, -otudf$sums), ]
    
    # Head of otudf to check prevalence & sums
    cat("Top OTUs based on prevalence and sums:\n")
    print(head(otudf))
    
    # Set a prevalence threshold dynamically based on the number of reads/min to max
    prevalence_threshold <- 0.1 * reads  # Example: threshold is 10% of reads
    nOTUs <- sum(otudf$prev >= prevalence_threshold)
    
    # Threshold and number of OTUs to keep
    cat("Prevalence threshold (minimum samples an OTU must appear in):", prevalence_threshold, "\n")
    cat("Number of OTUs to keep:", nOTUs, "\n")
    
    # Prune taxa based on prevalence
    physeq_16S_Prev_adj <- phyloseq::prune_taxa(rownames(otudf)[1:nOTUs], spiked_16S)
    
    # Save the adjusted prevalence phyloseq object
    saveRDS(physeq_16S_Prev_adj, file = output_file)
    cat("Adjusted prevalence phyloseq object saved to:", output_file, "\n")
    
    return(physeq_16S_Prev_adj)
  })
}

# Examples:
# Adjust prevalence based on the desired method
# spiked_16S_min <- adjusted_prevalence(physeq_16SOTU, method = "min")
# spiked_16S_max <- adjusted_prevalence(physeq_16SOTU, method = "max")
