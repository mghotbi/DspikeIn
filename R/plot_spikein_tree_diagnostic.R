#' @title Spike-in Tree Diagnostic Plot
#'
#' @description
#' Diagnostic visualization of spike-in taxa (ASVs/OTUs) on a phylogenetic tree.
#' Works with both `phyloseq` and `TreeSummarizedExperiment` objects (auto-converts internally).
#'
#' Visual elements:
#' - Tip labels (OTU/ASV names)
#' - Branch length annotations
#' - Prevalence (star size on tip)
#' - Log10(mean abundance) (bar ring)
#' - Sample metadata (discrete tile at tip)
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object filtered to spike-in taxa.
#' @param metadata_var Character. Sample metadata variable to use as tile color.
#' @param output_prefix Character. Filename prefix for saved plot. Default = "spikein_diag".
#' @param layout Character. Tree layout. Default = "circular".
#' @param save_plot Logical. Save figure if TRUE.
#' @param width, height Numeric. Width and height in inches.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @importFrom phyloseq psmelt phy_tree sample_data
#' @importFrom dplyr group_by summarise ungroup distinct
#' @importFrom ggplot2 aes theme ggsave scale_fill_manual scale_alpha_continuous scale_size_continuous scale_fill_gradient
#' @importFrom ggtree ggtree theme_tree2 geom_tiplab geom_text2
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom DspikeIn convert_tse_to_phyloseq
#' @importFrom ggstar geom_star
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load test data (16S OTU table with taxonomy, tree, and metadata)
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Filter to a known spike-in genus
#'   spikein <- phyloseq::subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
#'
#'   # Plot using phyloseq object (preview only)
#'   plot_spikein_tree_diagnostic(
#'     obj = spikein,
#'     metadata_var = "Animal.type",
#'     save_plot = FALSE
#'   )
#'
#'   # TreeSummarizedExperiment object
#'   # To save the plot
#'   # Use save_plot = TRUE & output_prefix = "tetragenococcus_diag"
#'
#'   tse_spikein <- convert_phyloseq_to_tse(spikein)
#'
#'   # Generate and save tree diagnostic plot
#'   plot_spikein_tree_diagnostic(
#'     obj = tse_spikein,
#'     metadata_var = "Animal.type",
#'     save_plot = FALSE,
#'   )
#' }
#' @export
plot_spikein_tree_diagnostic <- function(obj,
                                         metadata_var,
                                         output_prefix = "spikein_diag",
                                         layout = "circular",
                                         save_plot = FALSE,
                                         width = 10,
                                         height = 10) {
  # --- Accept both phyloseq and TSE objects ---
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  if (!inherits(obj, "phyloseq")) stop("Input must be a phyloseq object or convertible from TSE.")
  if (!metadata_var %in% colnames(sample_data(obj))) {
    stop("metadata_var not found in sample_data().")
  }

  # --- Prepare data ---
  df <- phyloseq::psmelt(obj)
  df$group <- df[[metadata_var]]

  df_summary <- df |>
    dplyr::group_by(OTU) |>
    dplyr::summarise(
      prevalence = mean(Abundance > 0),
      log_mean_abundance = log10(mean(Abundance) + 1),
      .groups = "drop"
    )

  # --- Base tree ---
  p <- ggtree::ggtree(phyloseq::phy_tree(obj), layout = layout) +
    ggtree::theme_tree2()

  # --- Tip labels ---
  p <- p + ggtree::geom_tiplab(size = 2, align = TRUE, linesize = 0.25, offset = -0.09)

  # --- Branch lengths ---
  p <- p + ggtree::geom_text2(
    aes(label = round(branch.length, 2)),
    hjust = -0.4,
    size = 2.2,
    color = "grey30"
  )

  # --- Prevalence (Star) ---
  p <- p + ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = df_summary,
      geom = geom_star,
      mapping = aes(y = OTU, x = 1.5, size = prevalence),
      starstroke = 0.1,
      fill = "#FB5607",
      pwidth = 0.05
    ) +
    ggplot2::scale_size_continuous(name = "Prevalence", range = c(1, 4))

  # --- log10(mean abundance) Bar ---
  p <- p + ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = df_summary,
      geom = geom_bar,
      mapping = aes(y = OTU, x = log_mean_abundance, fill = log_mean_abundance),
      orientation = "y",
      stat = "identity",
      pwidth = 0.1
    ) +
    ggplot2::scale_fill_gradient(
      name = "log10(mean abundance)",
      low = "white", high = "#3A86FF"
    )

  # --- Metadata tile ring ---
  p <- p + ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = dplyr::distinct(df, OTU, group),
      geom = geom_tile,
      mapping = aes(y = OTU, fill = group),
      width = 0.05,
      pwidth = 0.08
    ) +
    ggplot2::scale_fill_manual(
      values = DspikeIn::color_palette$mix_MG,
      name = metadata_var
    )

  # --- Save to file ---
  if (isTRUE(save_plot)) {
    ggplot2::ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
    ggplot2::ggsave(paste0(output_prefix, ".png"), p, width = width, height = height)
    message("Plots saved to: ", output_prefix, ".pdf/.png")
  }

  invisible(p)
}


# Example
# data("physeq_16SOTU", package = "DspikeIn")
#  spikein <- phyloseq::subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
#  spikein@sam_data$well.location

#  Sampale_spike<-plot_spikein_tree_diagnostic(
#   obj = spikein,
#   metadata_var = "Animal.type",
#   save_plot = TRUE,
#   output_prefix = "tetragenococcus_diag"  )
