#' @title Generate a Taxa Barplot with Relative or Absolute Abundance
#'
#' @description This function creates a bar plot of relative or absolute abundances of the top `n` taxa
#' (OTUs/ASVs) at a specified taxonomic rank. Non-top taxa are aggregated into an "Others" category,
#' which is always displayed at the top for clarity. The function supports normalization, faceting,
#' and customization of plot aesthetics. The plot is built using `ggplot2`, and the taxa are ordered
#' with "Others" appearing first.
#' @source Built on public functions from phyloseq and ggplot2 for data transformation and plotting.
#'
#' @param physeq A `phyloseq` or `TreeSummarizedExperiment` object containing microbiome data.
#' @param target_glom A character string specifying the taxonomic rank to plot (e.g., "Genus").
#' @param custom_tax_names A character vector specifying custom taxonomic names for the levels. Default is `NULL`.
#' @param normalize A logical value indicating whether to normalize sample counts to relative abundances. Default is `TRUE`.
#' @param treatment_variable A character string specifying the treatment variable for the x-axis (e.g., "Treatment"). Default is `"Treatment"`.
#' @param abundance_type A character string specifying whether to plot `"relative"` or `"absolute"` abundance. Default is `"relative"`.
#' @param x_angle A numeric value specifying the angle of x-axis text labels. Default is `25`.
#' @param fill_variable A character string specifying the variable to use for filling bar colors. Default is the same as `target_glom`.
#' @param facet_variable A character string specifying the variable to use for faceting the plot. Default is `NULL` (no faceting).
#' @param top_n_taxa A numeric value specifying the number of top taxa to include in the plot. Default is `20`.
#' @param palette A character vector of color codes or a function generating such a palette. Default is `color_palette$MG`.
#' @param legend_size A numeric value specifying the size of the legend text. Default is `11`.
#' @param legend_columns A numeric value specifying the number of columns for the legend. Default is `1`.
#' @param x_scale A character string specifying the x-axis scale in facets. Options are `"free"` or `"fixed"`. Default is `"free"`.
#' @param xlab A character string specifying the x-axis label. Default is `NULL` (no label).
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{barplot}{A `ggplot2` object representing the taxa barplot.}
#'   \item{taxa_data}{A `phyloseq` object containing the top taxa and aggregated "Others" taxa.}
#' }
#' @importFrom phyloseq tax_glom prune_taxa tax_table otu_table sample_data psmelt transform_sample_counts
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme scale_y_continuous element_text element_line
#' @importFrom ggplot2 guide_legend guides facet_grid element_rect theme_minimal
#' @importFrom rlang sym
#'
#' @examples
#' # Example 1: Relative abundance barplot for Genus
#' data("physeq_16SOTU", package = "DspikeIn")
#' bp_rel <- taxa_barplot(
#'   physeq = physeq_16SOTU,
#'   target_glom = "Genus",
#'   fill_variable = "Family",
#'   treatment_variable = "Diet",
#'   abundance_type = "relative",
#'   top_n_taxa = 20,
#'   legend_size = 10,
#'   x_scale = "free",
#'   legend_columns = 1,
#'   palette = DspikeIn::color_palette$MG
#' )
#' print(bp_rel$barplot)
#'
#' # Example 2: Absolute abundance barplot for Family with faceting
#' data("physeq_ITSOTU", package = "DspikeIn")
#' tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
#' bp_abs <- taxa_barplot(
#'   physeq = tse_ITSOTU,
#'   target_glom = "Genus",
#'   treatment_variable = "Animal.type",
#'   abundance_type = "absolute",
#'   fill_variable = "Family",
#'   facet_variable = "Diet",
#'   top_n_taxa = 10,
#'   x_scale = "fixed",
#'   xlab = NULL,
#'   legend_columns = 1,
#'   x_angle = 25,
#'   palette = DspikeIn::color_palette$cool_MG,
#'   legend_size = 12
#' )
#' print(bp_abs$barplot)
#'
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL,
                         normalize = TRUE, treatment_variable = "Treatment",
                         abundance_type = "relative", x_angle = 25, fill_variable = target_glom,
                         facet_variable = NULL, top_n_taxa = 20, palette = color_palette$MG,
                         legend_size = 11, legend_columns = 1, x_scale = "free", xlab = NULL) {
  #  Convert TSE to phyloseq if needed
  if (inherits(physeq, "TreeSummarizedExperiment")) {
    physeq <- convert_tse_to_phyloseq(physeq)
  }

  # Step 1: Taxonomic grouping at the specified rank
  glom <- phyloseq::tax_glom(physeq, taxrank = target_glom)

  # Step 2: Remove taxa with zero abundance across all samples
  glom <- phyloseq::prune_taxa(phyloseq::taxa_sums(glom) > 0, glom)

  # Step 3: Select top taxa based on abundance
  if (top_n_taxa > 0) {
    top_taxa <- names(
      sort(phyloseq::taxa_sums(glom), decreasing = TRUE)[seq_len(top_n_taxa)]
    )
  } else {
    top_taxa <- character(0)
  }

  # Step 4: Identify "Others" and combine
  otu_matrix <- as(phyloseq::otu_table(glom), "matrix")
  other_taxa <- setdiff(rownames(otu_matrix), top_taxa)

  if (length(other_taxa) > 0) {
    # Summing up "Others"
    others_abundance <- colSums(otu_matrix[other_taxa, , drop = FALSE])
    others_row <- matrix(others_abundance, nrow = 1, ncol = ncol(otu_matrix))
    rownames(others_row) <- "Others"

    # Combine top taxa and "Others"
    otu_combined <- rbind(otu_matrix[top_taxa, , drop = FALSE], others_row)

    # Rebuild taxonomy table with "Others"
    tax_top <- phyloseq::tax_table(glom)[top_taxa, , drop = FALSE]
    others_tax <- matrix(rep("Others", ncol(tax_top)), nrow = 1, dimnames = list("Others", colnames(tax_top)))
    tax_combined <- rbind(tax_top, others_tax)
  } else {
    otu_combined <- otu_matrix[top_taxa, , drop = FALSE]
    tax_combined <- phyloseq::tax_table(glom)[top_taxa, , drop = FALSE]
  }

  # Step 5: Ensure proper order ("Others" on top (alws), descending for others)
  row_order <- c("Others", setdiff(rownames(otu_combined), "Others"))
  otu_combined <- otu_combined[row_order, , drop = FALSE]
  tax_combined <- tax_combined[row_order, , drop = FALSE]

  # Step 6: Re-create phyloseq object
  new_physeq <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_combined, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_combined),
    phyloseq::sample_data(glom)
  )

  # Step 7: Normalize abundances if requested
  if (normalize && abundance_type == "relative") {
    new_physeq <- phyloseq::transform_sample_counts(new_physeq, function(x) x / sum(x))
  }

  # Step 8: Melt data for ggplot
  pm <- phyloseq::psmelt(new_physeq)
  pm <- pm[!is.na(pm$Abundance), ]

  # Ensure proper order for "fill_variable"
  pm[[fill_variable]] <- factor(
    pm[[fill_variable]],
    levels = c("Others", setdiff(unique(pm[[fill_variable]]), "Others"))
  )

  # Step 9: Build ggplot
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

  # Add faceting if specified
  if (!is.null(facet_variable)) {
    facet_var <- rlang::sym(facet_variable)
    facet_scales <- ifelse(x_scale == "free", "free_x", "fixed")
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(!!facet_var), scales = facet_scales)
  }

  # Add or suppress x-axis label
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab)
  } else {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }

  return(list(barplot = p, taxa_data = new_physeq))
}


# Example with fixed scales
# # ps is phyloseq object
#
# bp_fixed <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   abundance_type = "absolute",
#   facet_variable = "Treatment",
#   fill_variable = "Genus",
#   legend_size = 10,
#   top_n_taxa = 30,
#   xlab = NULL,
#   legend_columns = 1,
#   x_angle = 25,
#   legend_columns = 1,
#   x_scale = "fixed",
#   palette = DspikeIn::color_palette$light_MG
# )
# print(bp_fixed$barplot)
#
# # Example with free scales
# bp_free <- taxa_barplot(
#   physeq = tse_data,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   fill_variable = "Genus",
#   abundance_type = "relative",
#   facet_variable = "Treatment",
#   legend_size = 10,
#   top_n_taxa = 30,
#   legend_columns = 1,
#   x_scale = "free",
#   palette = DspikeIn::color_palette$extended_palette
# )
# print(bp_free$barplot)
#
#
# Example with free x-axis scale
# bp_free <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   fill_variable = "Family",
#   abundance_type = "absolute",
#   facet_variable = "Treatment",
#   x_scale = "free",
#   legend_size = 10,
#   top_n_taxa = 40,
#   legend_columns = 1,
#   palette = color_palette$MG)
# print(bp_free$barplot)

# bp_free <- taxa_barplot(
#   physeq= phy_M_M_ps,
#  target_glom = "Genus",
#  treatment_variable = "Diet",
# fill_variable = "Family",
#   abundance_type = "absolute",
#  facet_variable = "Habitat",
#  x_scale = "free",
#  legend_size = 10,
#  top_n_taxa = 40,
#  legend_columns = 1,
#  palette = color_palette$light_MG)
# print(bp_free$barplot)
