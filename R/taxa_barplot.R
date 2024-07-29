#' MG Color Palette
#'
#' This function returns a character vector of color palettes used in the package.
#'
#' @return A character vector of color hex codes.
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
    "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream",  "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue"
  )
}
#' Generate a Taxa Barplot
#'
#' This function creates a bar plot of the relative or absolute abundances of taxa at a specified taxonomic rank.
#' The top taxa are selected, and the plot can be customized with various options.
#' The function uses the `tax_glom` function from the `phyloseq` package for taxonomic grouping.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param target_glom A character string specifying the taxonomic rank to plot (e.g., "Genus").
#'        This specifies the taxonomic level at which the data should be aggregated and visualized.
#' @param custom_tax_names A character vector specifying custom taxonomic names for the levels. Default is NULL.
#'        This allows the user to provide custom names for the taxonomic levels in the taxonomy table of the phyloseq object.
#' @param normalize A logical indicating whether to normalize the sample counts to relative abundances. Default is TRUE.
#' @param treatment_variable A character string specifying the treatment variable to use for the x-axis. Default is "Treatment".
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "relative".
#' @param x_angle A numeric value specifying the angle of the x-axis text labels. Default is 20.
#' @param fill_variable A character string specifying the variable to use for fill in stacking taxa. Default is target_glom.
#' @param facet_variable A character string specifying the variable to use for faceting. Default is "Phylum".
#' @param top_n_taxa A numeric value specifying the number of top taxa to include in the plot. Default is 20.
#' @return A list containing the ggplot2 bar plot object (`barplot`) and the pruned phyloseq object (`taxa_data`) with the top taxa.
#' @examples
#' # Generate a taxa barplot for the Genus rank with relative abundance
#' bp <- taxa_barplot(physeqASV16, target_glom = "Genus", normalize = TRUE, 
#' treatment_variable = "animal.type", abundance_type = "relative")
#' print(bp$barplot)
#'
#' # Generate a taxa barplot for the Genus rank with absolute abundance
#' bp <- taxa_barplot(physeqASV16, target_glom = "Genus", normalize = FALSE, 
#' treatment_variable = "animal.type", abundance_type = "absolute")
#' print(bp$barplot)
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, 
                         normalize = TRUE, treatment_variable = "Treatment", 
                         abundance_type = "relative", x_angle = 20, fill_variable = target_glom, 
                         facet_variable = "Phylum", top_n_taxa = 20) {
  
  # Taxonomic grouping using tax_glom from phyloseq
  glom <- tax_glom(physeq, taxrank = target_glom)
  glom_1 <- prune_taxa(taxa_sums(glom) > 0, glom)
  glom_C <- tax_glom(glom_1, taxrank = target_glom)
  
  # Select top taxa
  top_taxa <- names(sort(taxa_sums(glom_C), decreasing = TRUE)[1:top_n_taxa])
  top_taxa_pruned <- prune_taxa(top_taxa, glom_C)
  top_v5 <- prune_taxa(taxa_sums(top_taxa_pruned) > 0, top_taxa_pruned)
  
  # Rename taxonomic levels
  if (!is.null(custom_tax_names)) {
    colnames(tax_table(top_v5)) <- custom_tax_names
  } else {
    colnames(tax_table(top_v5)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", target_glom)
  }
  
  # Normalize sample counts if requested
  if (normalize && abundance_type == "relative") {
    top_v5 <- transform_sample_counts(top_v5, function(OTU) OTU / sum(OTU))
  }
  
  # Prepare the data for ggplot
  pm <- psmelt(top_v5)
  
  # Create barplot
  if (abundance_type == "relative") {
    p <- ggplot(pm, aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format())
  } else {
    p <- ggplot(pm, aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      geom_bar(stat = "identity")
  }
  
  # Customize the plot
  p <- p + ggplot2::scale_fill_manual(values = MG()) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 20, keyheight = ggplot2::unit(0.6, "lines"), keywidth = ggplot2::unit(0.6, "lines"))) +
    ggplot2::labs(x = "", y = if (abundance_type == "relative") "Relative Abundance" else "Absolute Abundance") +
    ggplot2::facet_grid(cols = ggplot2::vars(.data[[facet_variable]]), scales = "free") +
    ggplot2::theme(
      strip.text.x = ggplot2::element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.text.y = ggplot2::element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.background = ggplot2::element_blank(),  # Remove the border of facets
      axis.text.x = ggplot2::element_text(angle = x_angle, vjust = 0.5, hjust = 1)
    )
  
  return(list(barplot = p, taxa_data = top_taxa_pruned))
}

# Example usage:
# Generate a taxa barplot for the Genus rank with absolute abundance
# bp_ab <- taxa_barplot(SP_Plethodon, target_glom = "Genus",
# treatment_variable = "Host.species", abundance_type = "absolute", x_angle = 90,
# fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
# print(bp_ab$barplot)

# Generate a taxa barplot for the Genus rank with relative abundance
#bp_rel <- taxa_barplot(Salamander_relative_NospikeSp, target_glom = "Genus", 
#treatment_variable = "Host.genus", abundance_type = "relative", x_angle = 90, 
#fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
#print(bp_rel$barplot)