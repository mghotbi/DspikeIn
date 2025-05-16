#' @title Taxa Bar Plot Without Aggregation (Relative or Absolute Abundance)
#'
#' @description
#' Create a bar plot of relative or absolute abundances of microbial taxa
#' **without using taxonomic glomming** (`tax_glom`). Works on raw or normalized data.
#'
#' @param physeq A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param tax_level A character string indicating the taxonomic rank (e.g., `"Genus"`).
#' @param normalize Logical; if `TRUE`, transforms to relative abundance. Default: `TRUE`.
#' @param treatment_variable Character; column name in sample metadata for x-axis grouping.
#' @param abundance_type Character; either `"relative"` or `"absolute"`. Default: `"relative"`.
#' @param x_angle Numeric; angle of x-axis tick labels. Default: `25`.
#' @param fill_variable Character; variable to fill bars by. Default: same as `tax_level`.
#' @param facet_variable Optional; column name for faceting. Default: `NULL`.
#' @param palette A named vector of colors to use for `fill_variable`.
#' @param legend_size Numeric; legend text size. Default: `11`.
#' @param legend_columns Integer; number of legend columns. Default: `1`.
#' @param x_scale Character; either `"free"` or `"fixed"` for facet x-axis scale. Default: `"free"`.
#' @param xlab Optional; override x-axis label. If `NULL`, x label is hidden.
#'
#' @return A `ggplot2` object containing a bar plot of taxa abundance.
#'
#' @importFrom phyloseq prune_taxa tax_table otu_table sample_data psmelt transform_sample_counts
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme scale_y_continuous element_text element_line
#' @importFrom ggplot2 guide_legend guides facet_wrap element_rect theme_minimal xlab
#' @importFrom rlang sym
#' @importFrom dplyr filter
#'
#' @examples
#' # Load required package
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'
#'   # Load example data
#'   data("physeq_ITSOTU", package = "DspikeIn")
#'
#'   # Subset: Eurycea salamanders from Blue Ridge, exclude unwanted genera
#'   Des <- physeq_ITSOTU |>
#'     phyloseq::subset_taxa(Genus != "Dekkera") |>
#'     phyloseq::subset_samples(Clade.Order == "Caudate") |>
#'     phyloseq::subset_samples(Host.genus == "Eurycea") |>
#'     phyloseq::subset_samples(Ecoregion.III == "Blue Ridge")
#'
#'   # Clean taxa: remove NAs or blanks in Phylum, filter low-abundance
#'   Des_filtered <- phyloseq::subset_taxa(Des, !is.na(Phylum) & Phylum != "")
#'   Des_ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(Des_filtered) > 99, Des_filtered)
#'
#'   # Plot taxa abundance with full control
#'   plotbar_abundance(
#'     physeq = Des_ps,
#'     normalize = TRUE,
#'     treatment_variable = "Diet",
#'     abundance_type = "absolute",
#'     x_angle = 0,
#'     fill_variable = "Phylum",
#'     palette = DspikeIn::color_palette$mix_MG,
#'     legend_size = 10,
#'     legend_columns = 1,
#'     x_scale = "free",
#'     xlab = NULL
#'   )
#' }
#'
#' @export
plotbar_abundance <- function(physeq, tax_level = "Genus", normalize = TRUE,
                              treatment_variable = "Host.taxon", abundance_type = "relative",
                              x_angle = 25, fill_variable = tax_level,facet_variable= NULL,
                              palette = DspikeIn::color_palette$mix_MG, legend_size = 11,
                              legend_columns = 1, x_scale = "free", xlab = NULL) {
  # Convert `TreeSummarizedExperiment` to `phyloseq` if needed
  if (inherits(physeq, "TreeSummarizedExperiment")) {
    physeq <- convert_tse_to_phyloseq(physeq)
  }

  # Ensure OTU and Tax Table Match
  common_taxa <- intersect(taxa_names(otu_table(physeq)), taxa_names(tax_table(physeq)))
  physeq <- phyloseq::prune_taxa(common_taxa, physeq)

  # Normalize if requested
  if (normalize && abundance_type == "relative") {
    physeq <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x))
  }

  # convert Phyloseq to Long Format for ggplot2
  pm <- phyloseq::psmelt(physeq)

  # Ensure Faceting Variable Exists
  # First: Check facet_variable is valid BEFORE converting to symbol
  if (!is.null(facet_variable)) {
    if (!(facet_variable %in% colnames(pm))) {
      stop(paste("Facet variable", facet_variable, "not found in sample metadata!"))
    }
    
    # Now it's safe to convert to symbol for ggplot
    facet_var <- rlang::sym(facet_variable)
    facet_scales <- ifelse(x_scale == "free", "free_x", "fixed")
    p <- p + ggplot2::facet_wrap(~ !!facet_var, scales = facet_scales)
  }

  
  # Build ggplot2 Barplot
  treatment_var <- rlang::sym(treatment_variable)
  fill_var <- rlang::sym(fill_variable)

  p <- ggplot2::ggplot(pm, ggplot2::aes(x = !!treatment_var, y = Abundance, fill = !!fill_var)) +
    ggplot2::geom_bar(stat = "identity", position = ifelse(abundance_type == "relative", "fill", "stack")) +
    ggplot2::scale_y_continuous(name = ifelse(abundance_type == "relative", "Relative Abundance", "Absolute Abundance")) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.key.size = ggplot2::unit(0.5, "cm"),
      axis.text.x = ggplot2::element_text(angle = x_angle, vjust = 0.5, hjust = 1),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.8)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_columns))

  # Add Faceting (If Specified)
  if (!is.null(facet_variable)) {
    facet_var <- rlang::sym(facet_variable)
    facet_scales <- ifelse(x_scale == "free", "free_x", "fixed")
    p <- p + ggplot2::facet_wrap(~ !!facet_var, scales = facet_scales)
  }

  # Add or Suppress X-Axis Label
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab)
  } else {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }

  return(p)
}

# RelSig16 <- subset_taxa(RelSig16, Family != "Chloroplast" & Order != "Chloroplast")
#
# barsigrelITS<-plotbar_abundance(sigrelITS, tax_level = "Genus",
#                                 treatment_variable = "Host.taxon",
#                                 x_angle = 6, abundance_type = "relative")
# RelSig16bar <-plotbar_abundance(rel_16S_sigd, tax_level = "Genus",
#                                 treatment_variable = "Host.taxon",
#                                 x_angle = 6, abundance_type = "absolute")

