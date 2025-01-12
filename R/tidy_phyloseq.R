#' Tidy a Phyloseq Object with steps below
#'
#' This function cleans and tidies a phyloseq object by performing the following steps:
#' - Fixes taxa names by removing any characters followed by '__' and any spaces after '__'
#' - Sets taxonomic ranks to standard names and requires 7 ranks
#' - Trims leading and trailing whitespace from taxa names
#' - Removes taxa with zero counts across all samples
#' - Removes taxa classified as "Chloroplast" at the Class level
#' - Removes taxa classified as "Mitochondria" at the Family level
#'
#' @param my_phyloseq A phyloseq object containing the taxonomic and abundance data.
#' @return A cleaned and tidied phyloseq object.
#' @importFrom phyloseq tax_table prune_taxa taxa_sums subset_taxa otu_table
#' @examples
#' \dontrun{
#' # Example usage:
#' spiked_16S <- tidy_phyloseq(spiked_16S)
#' }
#' @export
tidy_phyloseq <- function(my_phyloseq) {
  
  # Fix taxa names by removing any characters followed by '__' and any spaces after '__'
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    phyloseq::tax_table(my_phyloseq)[, col] <- gsub("[a-z]__\\s*", "", phyloseq::tax_table(my_phyloseq)[, col])
  }
  
  # Set taxonomic ranks to standard 7 ranks (if they exist)
  tax_table_cols <- colnames(phyloseq::tax_table(my_phyloseq))
  required_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (all(required_ranks %in% tax_table_cols)) {
    phyloseq::tax_table(my_phyloseq) <- phyloseq::tax_table(my_phyloseq)[, required_ranks]
  } else {
    warning("Some standard taxonomic ranks (Kingdom to Species) are missing. Proceeding with available ranks.")
  }
  
  # Trim leading and trailing whitespace from taxa names
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    phyloseq::tax_table(my_phyloseq)[, col] <- trimws(phyloseq::tax_table(my_phyloseq)[, col])
  }
  
  # Remove taxa with zero counts across all samples
  zero_row_indices <- rowSums(phyloseq::otu_table(my_phyloseq)) == 0
  my_phyloseq <- phyloseq::prune_taxa(!zero_row_indices, my_phyloseq)
  
  # Remove taxa classified as "Chloroplast" at the Class level
  if ("Class" %in% colnames(phyloseq::tax_table(my_phyloseq))) {
    my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Class != "Chloroplast")
  } else {
    warning("The taxonomic rank 'Class' is not present in the tax_table. Skipping removal of 'Chloroplast'.")
  }
  
  # Remove taxa classified as "Mitochondria" at the Family level
  if ("Family" %in% colnames(phyloseq::tax_table(my_phyloseq))) {
    my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Family != "Mitochondria")
  } else {
    warning("The taxonomic rank 'Family' is not present in the tax_table. Skipping removal of 'Mitochondria'.")
  }
  
  return(my_phyloseq)
}

# Example usage:
# ps <- tidy_phyloseq(ps)
