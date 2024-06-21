#' Tidy a Phyloseq Object
#'
#' This function cleans and tidies a phyloseq object by performing the following steps:
#' - Fixes taxa names by removing any characters followed by '__' and any spaces after '__'
#' - Sets taxonomic ranks to standard names
#' - Trims leading and trailing whitespace from taxa names
#' - Replaces NA in the Phylum column with "Unidentified"
#' - Removes taxa with zero counts
#' - Removes taxa classified as "Chloroplast" at the Class level
#' - Removes taxa classified as "Mitochondria" at the Family level
#'
#' @param my_phyloseq A phyloseq object containing the taxonomic and abundance data.
#' @return A cleaned and tidied phyloseq object.
#' @examples
#' # Example usage:
#' # spiked_16S <- tidy_phyloseq(spiked_16S)
#' @export
tidy_phyloseq <- function(my_phyloseq) {
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  
  library(phyloseq)
  
  # Fix taxa names by removing any characters followed by '__' and any spaces after '__'
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    phyloseq::tax_table(my_phyloseq)[, col] <- gsub("[a-z]__\\s*", "", phyloseq::tax_table(my_phyloseq)[, col]) 
  }
  
  # Set taxonomic ranks
  colnames(phyloseq::tax_table(my_phyloseq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # Fix taxa names
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    # Trim leading and trailing whitespace from taxa names
    phyloseq::tax_table(my_phyloseq)[, col] <- trimws(phyloseq::tax_table(my_phyloseq)[, col])
  }
  
  # Replace NA in Phylum with "Unidentified"
  phyloseq::tax_table(my_phyloseq)[is.na(phyloseq::tax_table(my_phyloseq)[, "Phylum"]), "Phylum"] <- "Unidentified"
  
  # Remove taxa with zero counts
  if (phyloseq::ntaxa(my_phyloseq) > 0) {
    my_phyloseq <- phyloseq::prune_taxa(taxa_sums(my_phyloseq) > 0, my_phyloseq)
  }
  
  # Remove taxa classified as "Chloroplast" at the Class level
  my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Class != "Chloroplast")
  
  # Remove taxa classified as "Mitochondria" at the Family level
  my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Family != "Mitochondria")
  
  return(my_phyloseq)
}

# Example usage:
# spiked_16S <- tidy_phyloseq(spiked_16S)
