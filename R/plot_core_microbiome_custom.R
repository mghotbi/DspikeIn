#' Plot Core Microbiome Prevalence Heatmap
#'
#' This function generates a prevalence heatmap of the core microbiome at a specified taxonomic rank.
#' It includes a detection threshold for prevalence, performs glomming, pruning, and proportion transformation,
#' and renames taxa with an ASV prefix. The function requires the `phyloseq`, `microbiomeutilities`, and `ggplot2` packages.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param taxrank A character string specifying the taxonomic rank to glom taxa. Default is "Genus".
#' @param select_taxa A character vector of taxa to select. Default is NULL, meaning no specific taxa are selected.
#' @param prevalences A numeric vector specifying the prevalence thresholds for plotting. Default is seq(0.01, 1, 0.01).
#' @param detections A numeric vector specifying the detection thresholds for plotting. Default is 10^seq(log10(1e-3), log10(0.3), length = 5).
#' @param min_prevalence A numeric value specifying the minimum prevalence threshold for core microbiome. Default is 0.3.
#' @param output_core_csv A character string specifying the path to save the core microbiome subset as a CSV file. Default is NULL, meaning no CSV file is saved.
#' @param output_core_rds A character string specifying the path to save the core microbiome subset as an RDS file. Default is NULL, meaning no RDS file is saved.
#' @return A ggplot2 object representing the core microbiome prevalence heatmap.
#' @examples
#' # Example usage:
#' # plot_core_microbiome_custom(physeq, taxrank = "Genus", select_taxa = NULL)
#' @export
plot_core_microbiome_custom <- function(physeq = NULL, 
                                        taxrank = "Genus", 
                                        select_taxa = NULL, 
                                        prevalences = seq(0.01, 1, 0.01), 
                                        detections = 10^seq(log10(1e-3), log10(0.3), length = 5),
                                        min_prevalence = 0.3,
                                        output_core_csv = NULL, 
                                        output_core_rds = NULL) {
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("microbiomeutilities", quietly = TRUE)) {
    stop("Package 'microbiomeutilities' is required but not installed.")
  }
  
  library(phyloseq)
  library(ggplot2)
  library(microbiomeutilities)
  library(RColorBrewer)
  
  # Check for required input
  if (is.null(physeq)) {
    stop("Error: 'physeq' argument is required.")
  }
  
  # Glom taxa at specified taxonomic rank
  glom_phy <- phyloseq::tax_glom(physeq, taxrank = taxrank)
  
  # Prune taxa with no abundance if select_taxa is specified
  if (!is.null(select_taxa)) {
    glom_phy <- prune_taxa(select_taxa, glom_phy)
  }
  
  # Exclude the species level
  glom_phy <- subset_taxa(glom_phy, select = -Species)
  
  # Transform counts to relative abundance
  glom_phy <- transform_sample_counts(glom_phy, function(x) 100 * x / sum(x))
  
  # Rename taxa with ASV prefix
  taxa_names(glom_phy) <- paste0("ASV", seq(ntaxa(glom_phy)))
  
  # Format phyloseq object for besthit
  phy_rel_f <- microbiomeutilities::format_to_besthit(glom_phy)
  
  # Save core microbiome as CSV file if output_core_csv argument is provided
  if (!is.null(output_core_csv)) {
    pm_core <- psmelt(phy_rel_f)
    write.csv(pm_core, output_core_csv, row.names = FALSE)
    cat("Core microbiome saved as CSV to:", output_core_csv, "\n")
  }
  
  # Save core microbiome as RDS file if output_core_rds argument is provided
  if (!is.null(output_core_rds)) {
    saveRDS(phy_rel_f, file = output_core_rds)
    cat("Core microbiome saved as RDS to:", output_core_rds, "\n")
  }
  
  # Plot core microbiome
  p.core <- plot_core(phy_rel_f, 
                      plot.type = "heatmap", 
                      colours = rev(brewer.pal(5, "Spectral")),
                      prevalences = prevalences, 
                      detections = detections, 
                      min.prevalence = min_prevalence) + 
    xlab("Detection Threshold (Relative Abundance (%))") +
    theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = 'white', color = "#e1deda"),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = 'black', size = 0.6),
      axis.line.y = element_line(colour = 'black', size = 0.6),
      axis.ticks = element_line(colour = 'black', size = 0.35),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12, color = "black", face = "bold"), 
      legend.key.size = unit(1, 'cm'),
      axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "plain"), 
      axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "plain"), 
      axis.text.x = element_text(family = "Times New Roman", size = 12, angle = 10, color = "black", face = "bold"), 
      axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      plot.title = element_text(color = "black", size = 12, face = "bold"),
      plot.subtitle = element_text(size = 11),
      strip.text.x = element_text(family="Times New Roman", size = 12, color="black", face = "bold"),
      strip.text.y = element_text(family="Times New Roman", size = 12, color="black", face = "bold"),
      strip.background = element_rect(colour = "black", fill = "gray98")
    )
  
  return(p.core)
}


# Example usage:
# Specify a different detection threshold
# custom_detections <- 10^seq(log10(3e-1), log10(0.5), length = 5)
# PCM <- plot_core_microbiome_custom(ps, detections = custom_detections, taxrank = "Family", output_core_rds = "core_microbiome.rds", output_core_csv = "core_microbiome.csv")
