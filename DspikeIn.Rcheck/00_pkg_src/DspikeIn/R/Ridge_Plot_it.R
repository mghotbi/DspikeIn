#' @title Generate Ridge Plots for Taxonomic Abundance Distribution
#' @description This function processes microbiome data (`phyloseq` or `TreeSummarizedExperiment`),
#' rarefies the dataset, performs proportion transformation, and generates ridge plots
#' to visualize the distribution of relative abundances at a specified taxonomic rank.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing taxonomic and abundance data.
#' @param taxrank A character string specifying the taxonomic rank for glomming and plotting. Default is "Genus".
#' @param rarefaction_depth A numeric value specifying the rarefaction depth. If NULL, it is set to 90\%
#'        of the minimum sample sum (for `phyloseq`). Default is NULL.
#' @param top_n An integer specifying the number of top taxa to include in the plot. Default is 10.
#' @return A `ggplot2` object representing the ridge plot of the distribution of relative abundances.
#'
#' @details
#' This function:
#' - Rarefies the dataset (if `phyloseq`) to normalize sample depth.
#' - Extracts abundance and taxonomic data using `get_otu_table()` and `get_tax_table()`.
#' - Aggregates abundance data at the specified taxonomic rank.
#' - Selects the top `n` taxa by total abundance.
#' - Generates a ridge plot to visualize abundance distributions.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load phyloseq object
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   ridge_physeq <- ridge_plot_it(physeq_16SOTU, taxrank = "Family", top_n = 10)
#'
#'
#'   # convert phyloseq object to TSE
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   ridge_tse <- ridge_plot_it(tse_16SOTU, taxrank = "Family", top_n = 10)
#' }
#'
#' @importFrom phyloseq rarefy_even_depth
#' @importFrom ggplot2 ggplot aes ggtitle labs theme_minimal theme element_text
#' @importFrom dplyr group_by summarise arrange desc slice_head pull filter
#' @importFrom rlang sym
#' @importFrom ggridges geom_density_ridges2
#' @importFrom reshape2 melt
#' @export
ridge_plot_it <- function(obj, taxrank = "Genus", rarefaction_depth = NULL, top_n = 10) {
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  if (inherits(obj, "phyloseq") && is.null(rarefaction_depth)) {
    min_depth <- min(phyloseq::sample_sums(obj))
    if (min_depth < 10) warning("Some samples have very low depth (<10); rarefaction may not be reliable.")
    rarefaction_depth <- max(10, floor(0.9 * min_depth))
    obj <- phyloseq::rarefy_even_depth(obj, rngseed = 500, sample.size = rarefaction_depth, replace = TRUE)
  }

  # Extract OTU and taxonomy tables using accessors
  otu_table <- get_otu_table(obj)
  tax_data <- get_tax_table(obj)

  # Ensure taxonomic rank exists
  if (is.null(tax_data) || !(taxrank %in% colnames(tax_data))) {
    stop("Taxonomic rank ", taxrank, " not found in dataset.")
  }

  # Aggregate abundance at the selected taxonomic level
  otu_table_aggregated <- aggregate(otu_table, by = list(tax_data[[taxrank]]), FUN = sum)
  colnames(otu_table_aggregated)[1] <- taxrank # Rename first column

  # Summarize total abundance per taxon and select top taxa
  top_taxa <- otu_table_aggregated |>
    dplyr::group_by(!!rlang::sym(taxrank)) |>
    dplyr::summarise(TotalAbundance = sum(dplyr::across(where(is.numeric)))) |>
    dplyr::arrange(dplyr::desc(TotalAbundance)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(!!rlang::sym(taxrank))

  # Filter dataset for selected top taxa
  otu_table_filtered <- otu_table_aggregated |>
    dplyr::filter(!!rlang::sym(taxrank) %in% top_taxa)

  # Prepare data for plotting
  plot_data <- reshape2::melt(otu_table_filtered, id.vars = taxrank, variable.name = "Sample", value.name = "Abundance")

  # Generate ridge plot
  plot <- ggplot2::ggplot(
    plot_data |> dplyr::filter(Abundance > 0 & !is.na(!!rlang::sym(taxrank))),
    ggplot2::aes(y = !!rlang::sym(taxrank), x = log10(Abundance), fill = !!rlang::sym(taxrank))
  ) +
    ggridges::geom_density_ridges2(scale = 1, alpha = 0.8, show.legend = FALSE) +
    ggplot2::ggtitle("Abundance Distribution") +
    ggplot2::labs(x = "Log10 (Abundance)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))

  return(plot)
}

# Usage Example
# ridge_tse <- ridge_plot_it(M19_tse, taxrank = "Family", top_n = 10)
# ridge_phy <- ridge_plot_it(M19, taxrank = "Family", top_n = 10)
