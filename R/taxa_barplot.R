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
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, normalize = TRUE, treatment_variable = "Treatment", abundance_type = "relative", x_angle = 20, fill_variable = target_glom, facet_variable = "Phylum", top_n_taxa = 20) {
  
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
  p <- p + scale_fill_manual(values = MG) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 20, keyheight = unit(0.6, "lines"), keywidth = unit(0.6, "lines"))) +
    my_custom_theme() +
    labs(x = "", y = if (abundance_type == "relative") "Relative Abundance" else "Absolute Abundance") +
    facet_grid(cols = vars(.data[[facet_variable]]), scales = "free") +
    theme(
      strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.background = element_blank(),  # Remove the border of facets
      axis.text.x = element_text(angle = x_angle, vjust = 0.5, hjust = 1)
    )
  
  return(list(barplot = p, taxa_data = top_v5))
}

# Example usage:
# Generate a taxa barplot for the Genus rank with absolute abundance
#bp_ab <- taxa_barplot(Salamander_absolute_NospikeSp, target_glom = "Genus",
#treatment_variable = "Host.genus", abundance_type = "absolute", x_angle = 90, 
#fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
#print(bp_ab$barplot)

# Generate a taxa barplot for the Genus rank with relative abundance
#bp_rel <- taxa_barplot(Salamander_relative_NospikeSp, target_glom = "Genus", 
#treatment_variable = "Host.genus", abundance_type = "relative", x_angle = 90, 
#fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
#print(bp_rel$barplot)