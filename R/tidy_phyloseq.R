#' Tidy a Phyloseq Object
#'
#' This function cleans and tidies a phyloseq object by performing the following steps:
#' - Fixes taxa names by removing any characters followed by '__' and any spaces after '__'
#' - Sets taxonomic ranks to standard names and requires 7 ranks
#' - Trims leading and trailing whitespace from taxa names
#' - Replaces NA in the Phylum column with "Unidentified"
#' - Removes taxa with zero counts
#' - Removes taxa classified as "Chloroplast" at the Class level
#' - Removes taxa classified as "Mitochondria" at the Family level
#'
#' @param my_phyloseq A phyloseq object containing the taxonomic and abundance data.
#' @return A cleaned and tidied phyloseq object.
#' @importFrom phyloseq tax_table prune_taxa taxa_sums subset_taxa
#' @examples
#' \dontrun{
#' # Note: The DspikeIn package functions only work with datasets 
#' # that contain 7 taxonomic ranks (e.g., Kingdom, Phylum, Class,
#' # Order, Family, Genus, Species).
#' 
#' # Example usage:
#' #' # Tidy up the phyloseq object for easier downstream analysis
#' spiked_16S <- tidy_phyloseq(spiked_16S)
#'
#' # The tidied phyloseq object can now be used in further analyses
#' }
#' @export
tidy_phyloseq <- function(my_phyloseq) {
  # Load necessary libraries with suppressMessages
  suppressMessages(if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  })
  
  # Fix taxa names by removing any characters followed by '__' and any spaces after '__'
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    phyloseq::tax_table(my_phyloseq)[, col] <- gsub("[a-z]__\\s*", "", phyloseq::tax_table(my_phyloseq)[, col])
  }
  
  # Set taxonomic ranks
  phyloseq::tax_table(my_phyloseq) <- phyloseq::tax_table(my_phyloseq)[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
  
  # Trim leading and trailing whitespace from taxa names
  for (col in colnames(phyloseq::tax_table(my_phyloseq))) {
    phyloseq::tax_table(my_phyloseq)[, col] <- trimws(phyloseq::tax_table(my_phyloseq)[, col])
  }
  
  # Replace NA in Phylum with "Unidentified"
  phyloseq::tax_table(my_phyloseq)[is.na(phyloseq::tax_table(my_phyloseq)[, "Phylum"]), "Phylum"] <- "Unidentified"
  
  # Remove taxa with zero counts
  my_phyloseq <- phyloseq::prune_taxa(phyloseq::taxa_sums(my_phyloseq) > 0, my_phyloseq)
  
  # Check if "Class" and "Family" columns exist before attempting to subset taxa
  if ("Class" %in% colnames(phyloseq::tax_table(my_phyloseq))) {
    # Remove taxa classified as "Chloroplast" at the Class level
    my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Class != "Chloroplast")
  } else {
    warning("The taxonomic rank 'Class' is not present in the tax_table. Skipping removal of 'Chloroplast'.")
  }
  
  if ("Family" %in% colnames(phyloseq::tax_table(my_phyloseq))) {
    # Remove taxa classified as "Mitochondria" at the Family level
    my_phyloseq <- phyloseq::subset_taxa(my_phyloseq, Family != "Mitochondria")
  } else {
    warning("The taxonomic rank 'Family' is not present in the tax_table. Skipping removal of 'Mitochondria'.")
  }
  
  return(my_phyloseq)
}

# Example usage:
# physeq_ITSOTU <- tidy_phyloseq(physeq_ITSOTU)
