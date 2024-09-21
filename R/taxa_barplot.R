#' Create a Taxa Barplot
#'
#' This function generates a bar plot of the relative or absolute abundances of taxa at a specified taxonomic rank.
#' The top taxa are selected, and the plot can be customized with various options.
#'
#' @importFrom phyloseq tax_glom prune_taxa tax_table taxa_sums psmelt transform_sample_counts
#' @importFrom ggplot2 ggplot geom_bar scale_y_continuous scale_fill_manual theme guides guide_legend labs facet_grid vars element_text element_blank element_line
#' @importFrom scales percent_format
#' @importFrom magrittr %>%
#' 
#' @param physeq A phyloseq object containing the microbiome data.
#' @param target_glom A character string specifying the taxonomic rank to plot (e.g., "Genus").
#' @param custom_tax_names A character vector specifying custom taxonomic names for the levels. Default is NULL.
#' @param normalize A logical indicating whether to normalize the sample counts to relative abundances. Default is TRUE.
#' @param treatment_variable A character string specifying the treatment variable to use for the x-axis. Default is "Treatment".
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "relative".
#' @param x_angle A numeric value specifying the angle of the x-axis text labels. Default is 20.
#' @param fill_variable A character string specifying the variable to use for fill in stacking taxa. Default is target_glom.
#' @param facet_variable A character string specifying the variable to use for faceting. Default is "Phylum".
#' @param top_n_taxa A numeric value specifying the number of top taxa to include in the plot. Default is 20.
#' @param palette A character vector of color hex codes to use for the fill colors. Default is `MG()`.
#'
#' @return A list containing the ggplot2 bar plot object (`barplot`) and the pruned phyloseq object (`taxa_data`) with the top taxa.
#'
#' @examples
#' \dontrun{
#' # Generate a taxa barplot for the Genus rank with absolute abundance
#' bp_ab <- taxa_barplot(physeq_16SASV, 
#'                       target_glom = "Genus",
#'                       treatment_variable = "Host.species", 
#'                       abundance_type = "absolute", 
#'                       x_angle = 90,
#'                       fill_variable = "Genus", 
#'                       facet_variable = "Diet", 
#'                       top_n_taxa = 20, 
#'                       palette = MG())
#' 
#' # Print the barplot for absolute abundance
#' print(bp_ab$barplot)
#'
#' # Generate a taxa barplot for the Genus rank with relative abundance
#' bp_rel <- taxa_barplot(physeq_16SASV, 
#'                        target_glom = "Genus", 
#'                        treatment_variable = "Host.genus", 
#'                        abundance_type = "relative", 
#'                        x_angle = 90, 
#'                        fill_variable = "Genus", 
#'                        facet_variable = "Diet", 
#'                        top_n_taxa = 20, 
#'                        palette = MG())
#'
#' # Print the barplot for relative abundance
#' print(bp_rel$barplot)
#' }
#' 
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, 
                         normalize = TRUE, treatment_variable = "Treatment", 
                         abundance_type = "relative", x_angle = 25, fill_variable = target_glom, 
                         facet_variable = "Phylum", top_n_taxa = 20, palette = MG()) {
  
  # Suppress specific package startup messages
  suppressWarnings({
    suppressPackageStartupMessages({
      if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop("Package 'phyloseq' is required but not installed.")
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
      }
    })
  })
  
  # Taxonomic grouping using tax_glom from phyloseq
  glom <- phyloseq::tax_glom(physeq, taxrank = target_glom)
  glom_1 <- phyloseq::prune_taxa(phyloseq::taxa_sums(glom) > 0, glom)
  glom_C <- phyloseq::tax_glom(glom_1, taxrank = target_glom)
  
  # Select top taxa
  top_taxa <- names(sort(phyloseq::taxa_sums(glom_C), decreasing = TRUE)[1:top_n_taxa])
  top_taxa_pruned <- phyloseq::prune_taxa(top_taxa, glom_C)
  top_v5 <- phyloseq::prune_taxa(phyloseq::taxa_sums(top_taxa_pruned) > 0, top_taxa_pruned)
  
  # Rename taxonomic levels
  if (!is.null(custom_tax_names)) {
    colnames(phyloseq::tax_table(top_v5)) <- custom_tax_names
  } else {
    colnames(phyloseq::tax_table(top_v5)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", target_glom)
  }
  
  # Normalize sample counts if requested
  if (normalize && abundance_type == "relative") {
    top_v5 <- phyloseq::transform_sample_counts(top_v5, function(OTU) OTU / sum(OTU))
  }
  
  # Prepare the data for ggplot
  pm <- phyloseq::psmelt(top_v5)
  
  # Create barplot
  if (abundance_type == "relative") {
    p <- ggplot2::ggplot(pm, ggplot2::aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      ggplot2::geom_bar(stat = "identity", position = "fill") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(), 
                                  name = "Relative Abundance (%)")  # Change Y-axis label
  } else {
    p <- ggplot2::ggplot(pm, ggplot2::aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::ylab("Absolute Abundance")
  }
  
  # Customize the plot with the specified palette and bold font, and add axis lines
  p <- p + 
    ggplot2::scale_fill_manual(values = palette) + # Custom color palette
    ggplot2::theme_minimal(base_size = 14) +  # Cleaner theme with larger text
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.height = ggplot2::unit(0.6, "lines"), 
      legend.key.width = ggplot2::unit(0.6, "lines"),
      axis.text.x = ggplot2::element_text(angle = x_angle, vjust = 0.5, hjust = 1),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.x = ggplot2::element_blank(),  # Removing x-axis title for cleaner look
      axis.text = ggplot2::element_text(size = 12),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.8),  # Add x-axis line
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.8),  # Add y-axis line
      panel.grid.major.x = ggplot2::element_blank(),  # Remove vertical gridlines
      panel.grid.minor = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(family = "Arial", size = 12, color = "black", face = "bold"),  # Ensure bold font
      strip.text.y = ggplot2::element_text(family = "Arial", size = 12, color = "black", face = "bold"),  # Ensure bold font
      strip.background = ggplot2::element_rect(fill = "gray90", color = NA)  # Light background for facet labels
    ) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data[[facet_variable]]), scales = "free")
  
  return(list(barplot = p, taxa_data = top_taxa_pruned))
}

#' MG Color Palette
#'
#' This function returns a character vector of color palettes used in the package.
#'
#' @return A character vector of color hex codes.
#' @examples
#' # Example usage of the MG color palette
#' my_colors <- MG()
#' print(my_colors)
#' @export
MG <- function() {
  c(
    "#FFFF33", "#FF7F00", "#E41A1C", "firebrick4", "#2e4057", "#984EA3", "#377EB8", "olivedrab3",
    "#4DAF4A", "#336633", "grey80", "#BB650B", "gold", "#559999", "#7570b3", "#E78AC3", "#A6D854",
    "#66a61e", "#e6ab02", "#a6761d", "#663300", "#66C2A5", "#0e669b", "#00798c", "dodgerblue4",
    "steelblue2", "#00AFBB", "#E7B800", "#FC4E07", "lightskyblue4", "green", "red", "#FFF000",
    "#0099CC", "#FF9933", "#CC9900", "chartreuse1", "#FF3399", "#00FFFF", "#0000CC", "#A37F6F",
    "#9183E6", "#00AD9A", "#990033", "#909800", "#00FF00", "#17b5b4", "#AED1D6", "#b1010c",
    "firebrick2", "blue", "navy", "yellow", "brown", "black", "purple", "darkred", "darkgreen",
    "#82cfd0", "#b2e0e4", "honeydew3", "#8d96a3", "lavender", "#CC6686", "lavenderblush2",
    "mistyrose3", "#e1deda", "darkgoldenrod", "burlywood", "papayawhip", "wheat4", "cornsilk3",
    "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream", 
    "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue"
  )
}

# Example usage:
# Generate a taxa barplot for the Genus rank with relative abundance
# In this example, we're grouping by "Host.genus" and faceting by "Diet"
# bp_rel <- taxa_barplot(physeq_16SASV, 
#                        target_glom = "Genus",            # Taxonomic level to glom
#                        treatment_variable = "Host.genus", # Variable to use on x-axis
#                        abundance_type = "relative",     # Plot relative 
#                        x_angle = 90,                   # Rotate x-axis 
#                        fill_variable = "Genus",       # Fill bars by Genus
#                        facet_variable = "Diet",      # Facet the plot by Diet
#                        top_n_taxa = 20,              # Show the top 20 taxa
#                        palette = MG())               # Use custom color palette (MG)
#
# print(bp_rel$barplot)
#
# Generate a taxa barplot for the Genus rank with absolute abundance
# In this example, we're grouping by "Host.species" and faceting by "Diet"
# bp_ab <- taxa_barplot(physeq_16SASV, 
#                       target_glom = "Genus",              # Taxonomic level to glom
#                       treatment_variable = "Host.species", # Variable to use on x-axis
#                       abundance_type = "absolute",    # Plot absolute abundance
#                       x_angle = 90,                   # Rotate x-axis labels by 90 degrees
#                       fill_variable = "Genus",         # Fill bars by Genus
#                       facet_variable = "Diet",       # Facet the plot by Diet
#                       top_n_taxa = 20,              # Show the top 20 taxa
#                       palette = MG())              # Use custom color palette (MG)
#
# print(bp_ab$barplot)
