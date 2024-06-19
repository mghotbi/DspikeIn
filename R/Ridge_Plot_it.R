#' Generate Ridge Plots to Visualize the Distribution of Abundances of Taxa
#'
#' This function rarefies the data, performs proportion transformation, gloms and prunes taxa, and generates ridge plots to visualize the distribution of relative abundances of taxa.
#'
#' @param physeq A phyloseq object containing the taxonomic and abundance data.
#' @param taxrank A character string specifying the taxonomic rank to use for glomming and plotting. Default is "Genus".
#' @param rarefaction_depth A numeric value specifying the rarefaction depth. If NULL, it is set to 90% of the minimum sample sum. Default is NULL.
#' @param top_n An integer specifying the number of top taxa to include in the plot. Default is 10.
#' @return A ggplot object representing the ridge plot of the distribution of relative abundances of taxa.
#' @examples
#' # Example usage:
#' # ridgeP <- ridge_plot_it(spiked_16S, taxrank = "Family", top_n = 10)
#' # print(ridgeP) + my_custom_theme()
#' @export
ridge_plot_it <- function(physeq, taxrank = "Genus", rarefaction_depth = NULL, top_n = 10) {
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  
  library(phyloseq)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  
  # Rarefy and tidy phyloseq object
  if (is.null(rarefaction_depth)) {
    rarefaction_depth <- max(10, 0.9 * min(sample_sums(physeq)))
  }
  physeq <- rarefy_even_depth(physeq, rngseed = 500, sample.size = rarefaction_depth, replace = TRUE)
  physeq <- tidy_phyloseq(physeq)
  
  # Glom and prune taxa, and filter out taxa with very small abundances
  glom <- tax_glom(physeq, taxrank = taxrank)
  glom <- prune_taxa(taxa_sums(glom) > 0, glom)  # Filter out taxa with very small abundances
  
  # Convert to data frame
  ps <- psmelt(glom)
  
  # Summarize and filter the top N taxa by abundance
  top_taxa <- ps %>%
    group_by_at(taxrank) %>%
    summarise(TotalAbundance = sum(Abundance)) %>%
    arrange(desc(TotalAbundance)) %>%
    head(top_n) %>%
    pull(!!sym(taxrank))
  
  ps <- ps %>% filter(!!sym(taxrank) %in% top_taxa)
  
  # Plot
  plot <- ggplot(ps %>% filter(Abundance > 0 & !is.na(get(taxrank))), aes_string(y = taxrank, x = "log10(Abundance)", fill = taxrank)) +
    geom_density_ridges2(scale = 1, alpha = 0.8, show.legend = FALSE) +
    ggtitle("Distribution of Relative Abundances") +
    labs(x = "Log10 (Abundance)", y = NULL) +
    theme_minimal() +
    scale_fill_manual(values = MG) +
    theme(axis.text.y = element_text(face = "italic"))
  
  return(plot)
}

# Example usage:
# ridgeP <- ridge_plot_it(spiked_16S, taxrank = "Family", top_n = 10)
# print(ridgeP) + my_custom_theme()

