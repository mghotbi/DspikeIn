#' Generate a Taxa Barplot with Relative or Absolute Abundance
#'
#' This function creates a bar plot of the relative or absolute abundances of the top `n` taxa 
#' (OTUs/ASVs) at a specified taxonomic rank. If desired, the function can aggregate 
#' non-top taxa into an "Others" category to ensure the bar fills up to 100% in relative abundance plots.
#'
#' @param physeq A phyloseq object containing microbiome data.
#' @param target_glom A character string specifying the taxonomic rank to plot (e.g., "Genus").
#' @param custom_tax_names A character vector specifying custom taxonomic names for the levels. Default is NULL.
#' @param normalize A logical value indicating whether to normalize the sample counts to relative abundances. Default is TRUE.
#' @param treatment_variable A character string specifying the treatment variable to use for the x-axis. Default is "Treatment".
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "relative".
#' @param x_angle A numeric value specifying the angle of the x-axis text labels. Default is 25.
#' @param fill_variable A character string specifying the variable to use for filling the bar colors. Default is the same as `target_glom`.
#' @param facet_variable A character string specifying the variable to use for faceting the plot. Default is NULL (no faceting).
#' @param top_n_taxa A numeric value specifying the number of top taxa to include in the plot. Default is 20.
#' @param palette A character vector of color codes or a function generating such a palette. Default is `color_palette$MG`.
#' @param legend_size A numeric value specifying the size of the legend text. Default is 11.
#' @param legend_columns A numeric value specifying the number of columns for the legend. Default is 1.
#' @param x_scale A character string specifying the x-axis scale in facets. Options are `"free"` or `"fixed"`. Default is `"fixed"`.
#' @param xlab A character string specifying the x-axis label. Default is `NULL` (no label).
#' 
#' @return A list containing the following components:
#' \describe{
#'   \item{barplot}{A ggplot2 object representing the taxa barplot.}
#'   \item{taxa_data}{A phyloseq object containing only the top taxa (and optionally "Others").}
#' }
#' 
#' @importFrom phyloseq tax_glom prune_taxa tax_table otu_table sample_data psmelt transform_sample_counts
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme scale_y_continuous element_text element_line
#' @importFrom ggplot2 guide_legend guides facet_grid element_rect theme_minimal
#' @importFrom rlang sym
#' 
#' @examples
#' \dontrun{
#' # Generate a barplot for the Genus rank with relative abundance
#' bp_rel <- taxa_barplot(
#'   physeq = physeq_data,                 # Your phyloseq object
#'   target_glom = "Genus",                # Taxonomic level to aggregate (e.g., Genus)
#'   treatment_variable = "Host.species",  # Variable to use on x-axis
#'   abundance_type = "relative",          # Plot relative abundance
#'   x_angle = 45,                         # Angle for x-axis labels
#'   fill_variable = "Genus",              # Fill bars by taxon (Genus)
#'   facet_variable = "Diet",              # Facet the plot by diet or NULL
#'   top_n_taxa = 20,                      # Show the top 20 taxa
#'   palette = color_palette$MG,           # Use the predefined MG color palette
#'   legend_size = 12,                     # Customize legend text size
#'   legend_columns = 2,                   # Split legend into 2 columns
#'   x_scale = "free",                     # Free x-axis levels across facets
#'   xlab = "Sample Groups"                # Add a custom x-axis label
#' )
#' print(bp_rel$barplot)                   # Print the resulting barplot
#' }
#' 
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, 
                         normalize = TRUE, treatment_variable = "Treatment", 
                         abundance_type = "relative", x_angle = 25, fill_variable = target_glom, 
                         facet_variable = NULL, top_n_taxa = 20, palette = color_palette$MG,
                         legend_size = 11, legend_columns = 1, x_scale = "fixed", xlab = NULL) {
  
  # Taxonomic grouping using tax_glom from phyloseq
  glom <- phyloseq::tax_glom(physeq, taxrank = target_glom)
  
  # Remove taxa with zero abundance across all samples
  glom_1 <- phyloseq::prune_taxa(phyloseq::taxa_sums(glom) > 0, glom)
  
  # Re-glom at the specified taxonomic level
  glom_C <- phyloseq::tax_glom(glom_1, taxrank = target_glom)
  
  # Select top taxa
  top_taxa <- names(sort(phyloseq::taxa_sums(glom_C), decreasing = TRUE)[1:top_n_taxa])
  top_taxa <- top_taxa[!is.na(top_taxa)]
  
  # Combine non-top taxa into "Others"
  otu_matrix <- as(phyloseq::otu_table(glom_C), "matrix")
  other_taxa <- setdiff(rownames(otu_matrix), top_taxa)
  if (length(other_taxa) > 0) {
    others_abundance <- colSums(otu_matrix[other_taxa, , drop = FALSE])
    others_row <- matrix(others_abundance, nrow = 1, ncol = ncol(otu_matrix))
    rownames(others_row) <- "Others"
    otu_top <- otu_matrix[top_taxa, , drop = FALSE]
    otu_combined <- rbind(otu_top, others_row)
    
    tax_top <- as(phyloseq::tax_table(glom_C)[top_taxa, ], "matrix")
    others_tax <- matrix(rep("Others", ncol(tax_top)), nrow = 1, dimnames = list("Others", colnames(tax_top)))
    tax_combined <- rbind(tax_top, others_tax)
  } else {
    otu_combined <- otu_matrix[top_taxa, , drop = FALSE]
    tax_combined <- as(phyloseq::tax_table(glom_C)[top_taxa, ], "matrix")
  }
  
  # Create new phyloseq object
  new_otu_table <- phyloseq::otu_table(otu_combined, taxa_are_rows = TRUE)
  new_tax_table <- phyloseq::tax_table(tax_combined)
  physeq_top_others <- phyloseq::phyloseq(new_otu_table, new_tax_table, phyloseq::sample_data(glom_C))
  
  # Rename taxonomic levels if custom names are provided
  if (!is.null(custom_tax_names)) {
    colnames(phyloseq::tax_table(physeq_top_others)) <- custom_tax_names
  } else {
    colnames(phyloseq::tax_table(physeq_top_others)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", target_glom)
  }
  
  # Normalize sample counts if requested
  if (normalize && abundance_type == "relative") {
    physeq_top_others <- phyloseq::transform_sample_counts(physeq_top_others, function(OTU) OTU / sum(OTU))
  }
  
  # Prepare the data for ggplot
  pm <- phyloseq::psmelt(physeq_top_others)
  treatment_var <- rlang::sym(treatment_variable)
  fill_var <- rlang::sym(fill_variable)
  
  # Create barplot
  p <- ggplot2::ggplot(pm, ggplot2::aes(x = !!treatment_var, y = Abundance, fill = !!fill_var)) +
    ggplot2::geom_bar(stat = "identity", position = ifelse(abundance_type == "relative", "fill", "stack")) +
    ggplot2::scale_y_continuous(name = ifelse(abundance_type == "relative", "Relative Abundance", "Absolute Abundance")) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.key.size = ggplot2::unit(0.5, "cm"), # Smaller legend symbols
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
  
  # Add or remove x-axis label
  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab)
  } else {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  
  return(list(barplot = p, taxa_data = physeq_top_others))
}



# # Example with fixed x-axis scale
# ps is phyloseq object
# bp_fix <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Family",
#   treatment_variable = "Genotype",
#   abundance_type = "absolute",
#   fill_variable = "Family",
#   facet_variable = "Treatment",
#   x_scale = "free",
#   legend_size = 10, 
#   top_n_taxa = 20, 
#   legend_columns = 1,
#   palette = color_palette$MG)
# print(bp_fix$barplot)
# 
# 
# bp_free <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   abundance_type = "absolute",
#   facet_variable = "Treatment",
#   fill_variable = "Genus",
#   x_scale = "free",
#   legend_size = 10, 
#   top_n_taxa = 30, 
#   legend_columns = 1,
#   xlab = NULL,
#   palette = color_palette$MG)
# print(bp_free$barplot)
# 
# 
# 
# bp_fix <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   abundance_type = "relative",
#   facet_variable = "Treatment",
#   fill_variable = "Genus",
#   x_scale = "fixed",
#   legend_size = 10,
#   top_n_taxa = 30,
#   xlab=NULL,
#   legend_columns = 2,
#   palette = color_palette$MG)
# print(bp_fix$barplot)

# 
# Example with free x-axis scale
# bp_free <- taxa_barplot(
#   physeq = ps,
#   target_glom = "Genus",
#   treatment_variable = "Genotype",
#   fill_variable = "Family",
#   abundance_type = "relative",
#   facet_variable = "Treatment",
#   x_scale = "free",
#   legend_size = 10,
#   top_n_taxa = 40,
#   legend_columns = 1,
#   palette = color_palette$MG)
# print(bp_free$barplot)
