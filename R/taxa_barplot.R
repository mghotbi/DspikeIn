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
#' @param palette Either a function that generates a color palette or a vector of color hex codes. Default is `MG()` from the package.
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
#' @importFrom dplyr %>%
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
#'   palette = color_palette$MG                       # Use a predefined color palette
#' )
#' print(bp_rel$barplot)                   # Print the resulting barplot
#' }
#' 
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, 
                         normalize = TRUE, treatment_variable = "Treatment", 
                         abundance_type = "relative", x_angle = 25, fill_variable = target_glom, 
                         facet_variable = NULL, top_n_taxa = 20, palette = MG()) {
  
  # Taxonomic grouping using tax_glom from phyloseq
  glom <- phyloseq::tax_glom(physeq, taxrank = target_glom)
  
  # Remove taxa with zero abundance across all samples
  glom_1 <- phyloseq::prune_taxa(phyloseq::taxa_sums(glom) > 0, glom)
  
  # Re-glom at the specified taxonomic level
  glom_C <- phyloseq::tax_glom(glom_1, taxrank = target_glom)
  
  # Select top taxa
  top_taxa <- names(sort(phyloseq::taxa_sums(glom_C), decreasing = TRUE)[1:top_n_taxa])
  
  # Remove any NA values from the top_taxa
  top_taxa <- top_taxa[!is.na(top_taxa)]
  
  # Extract OTU table and check if top_taxa exists in it
  otu_matrix <- as(phyloseq::otu_table(glom_C), "matrix")
  
  # Ensure that only top_taxa present in OTU matrix are selected
  top_taxa <- top_taxa[top_taxa %in% rownames(otu_matrix)]
  
  # Summing all non-top taxa as "Others" (only if there are remaining taxa)
  other_taxa <- setdiff(rownames(otu_matrix), top_taxa)
  
  # Conditionally add "Others" if there are any non-top taxa
  if (length(other_taxa) > 0) {
    others_abundance <- colSums(otu_matrix[other_taxa, , drop = FALSE])
    
    # Create a new "Others" taxon and merge it with top taxa
    others_row <- matrix(others_abundance, nrow = 1, ncol = ncol(otu_matrix))
    rownames(others_row) <- "Others"
    otu_top <- otu_matrix[top_taxa, , drop = FALSE]
    otu_combined <- rbind(otu_top, others_row)
    
    # Update the tax_table to include "Others"
    tax_top <- as(phyloseq::tax_table(glom_C)[top_taxa, ], "matrix")
    others_tax <- matrix(rep("Others", ncol(tax_top)), nrow = 1, dimnames = list("Others", colnames(tax_top)))
    tax_combined <- rbind(tax_top, others_tax)
  } else {
    # If no "Others", just use the top taxa
    otu_combined <- otu_matrix[top_taxa, , drop = FALSE]
    tax_combined <- as(phyloseq::tax_table(glom_C)[top_taxa, ], "matrix")
  }
  
  # Create a new phyloseq object with top taxa + "Others" (if applicable)
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
  
  # Convert dynamic variable names to symbols
  treatment_var <- rlang::sym(treatment_variable)
  fill_var <- rlang::sym(fill_variable)
  
  # Create barplot
  if (abundance_type == "relative") {
    p <- ggplot2::ggplot(pm, ggplot2::aes(x = !!treatment_var, y = Abundance, fill = !!fill_var)) +
      ggplot2::geom_bar(stat = "identity", position = "fill") +
      ggplot2::scale_y_continuous(name = "Relative Abundance")
  } else {
    p <- ggplot2::ggplot(pm, ggplot2::aes(x = !!treatment_var, y = Abundance, fill = !!fill_var)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::ylab("Absolute Abundance")
  }
  
  # Use either a color palette function or a color vector
  colors <- if (is.function(palette)) palette() else palette
  
  # Customize the plot with the specified palette and bold font
  p <- p + 
    ggplot2::scale_fill_manual(values = colors) +  # Use custom color palette
    ggplot2::theme_minimal(base_size = 14) +        # Cleaner theme with larger text
    ggplot2::theme(
      legend.position = "right",             
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = x_angle, vjust = 0.5, hjust = 1),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.8), 
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.8)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))  # Set number of columns to 1
  
  # Conditionally add faceting if facet_variable is provided
  if (!is.null(facet_variable)) {
    facet_var <- rlang::sym(facet_variable)
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(!!facet_var), scales = "free")
  }
  
  return(list(barplot = p, taxa_data = physeq_top_others))
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
#                        top_n_taxa = 20,              # Show the top 20 taxa
#                        palette = MG())               # Use custom color palette (MG)
#
# print(bp_rel$barplot)
#
# Generate a taxa barplot for the Genus rank with absolute abundance
# In this example, we're grouping by "Host.species" and faceting by "Diet"
# bp_ab <- taxa_barplot(ps, 
#                       target_glom = "Genus",              # Taxonomic level to glom
#                       treatment_variable = "Host.species", # Variable to use on x-axis
#                       abundance_type = "absolute",         # Plot absolute abundance
#                       x_angle = 90,                       # Rotate x-axis labels by 90 degrees
#                       fill_variable = "Genus",            # Fill bars by Genus
#                       facet_variable = "Diet",            # Facet the plot by Diet
#                       top_n_taxa = 20,                    # Show the top 20 taxa
#                       palette = color_palette$MG)                     # Use custom color palette (MG)
#
# print(bp_ab$barplot)                                      # Print the barplot
