#' Generate Ridge Plots to Visualize the Distribution of Abundances of Taxa
#'
#' This function rarefies the data, performs proportion transformation, gloms and prunes taxa, 
#' and generates ridge plots to visualize the distribution of relative abundances of taxa.
#'
#' @param physeq A phyloseq object containing the taxonomic and abundance data.
#' @param taxrank A character string specifying the taxonomic rank to use for glomming and plotting. Default is "Genus".
#' @param rarefaction_depth A numeric value specifying the rarefaction depth. If NULL, it is set to 90\% of the minimum sample sum. Default is NULL.
#' @param top_n An integer specifying the number of top taxa to include in the plot. Default is 10.
#' @return A ggplot2 object representing the ridge plot of the distribution of relative abundances of taxa.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ridgeP <- ridge_plot_it(spiked_16S, taxrank = "Family", top_n = 10)
#'   print(ridgeP) + my_custom_theme()
#' }
#' }
#' @importFrom phyloseq rarefy_even_depth tax_glom prune_taxa taxa_sums psmelt
#' @importFrom ggplot2 ggplot aes ggtitle labs theme_minimal theme element_text
#' @importFrom dplyr group_by summarise arrange desc slice_head pull filter
#' @importFrom rlang sym
#' @importFrom ggridges geom_density_ridges2
#' @export
ridge_plot_it <- function(physeq, taxrank = "Genus", rarefaction_depth = NULL, top_n = 10) {
  # Rarefy and tidy phyloseq object
  if (is.null(rarefaction_depth)) {
    rarefaction_depth <- max(10, 0.9 * min(phyloseq::sample_sums(physeq)))
  }
  physeq <- phyloseq::rarefy_even_depth(physeq, rngseed = 500, sample.size = rarefaction_depth, replace = TRUE)
  
  # Glom and prune taxa, and filter out taxa with very small abundances
  glom <- phyloseq::tax_glom(physeq, taxrank = taxrank)
  glom <- phyloseq::prune_taxa(phyloseq::taxa_sums(glom) > 0, glom)
  
  # Convert to data frame
  ps <- phyloseq::psmelt(glom)
  
  # Summarize and filter the top N taxa by abundance
  top_taxa <- ps %>%
    dplyr::group_by(!!rlang::sym(taxrank)) %>%
    dplyr::summarise(TotalAbundance = sum(Abundance)) %>%
    dplyr::arrange(dplyr::desc(TotalAbundance)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(!!rlang::sym(taxrank))
  
  ps <- ps %>% dplyr::filter(!!rlang::sym(taxrank) %in% top_taxa)
  
  # Plot
  plot <- ggplot2::ggplot(ps %>% dplyr::filter(Abundance > 0 & !is.na(!!rlang::sym(taxrank))),
                          ggplot2::aes(y = !!rlang::sym(taxrank), x = log10(Abundance), fill = !!rlang::sym(taxrank))) +
    ggridges::geom_density_ridges2(scale = 1, alpha = 0.8, show.legend = FALSE) +
    ggplot2::ggtitle("Abundance Distribution") +
    ggplot2::labs(x = "Log10 (Abundance)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
  
  return(plot)
}

# Example usage:
# ridgeP <- ridge_plot_it(physeq_16SOTU, taxrank = "Family", top_n = 10)
# print(ridgeP) + my_custom_theme()
