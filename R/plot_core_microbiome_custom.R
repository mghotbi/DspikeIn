#' Plot Core Microbiome Prevalence Heatmap
#'
#' This function generates a prevalence heatmap of the core microbiome at a specified taxonomic rank.
#' It allows users to pass custom detection thresholds, prevalence thresholds, and a minimum prevalence 
#' filter, and it provides an option to order taxa either in ascending or descending abundance.
#' The plot displays a heatmap showing detection thresholds at different prevalence levels.
#'
#' @param physeq A \code{phyloseq} object containing the microbial data.
#' @param taxrank A character string specifying the taxonomic rank to glom taxa. Default is "Genus".
#' @param select_taxa A character vector of taxa to select. Default is \code{NULL}, meaning no specific taxa are selected.
#' @param detections A list with the following elements:
#'   \itemize{
#'     \item \code{prevalences}: A numeric vector specifying the prevalence thresholds for plotting. Default is \code{seq(0.03, 1, 0.01)}.
#'     \item \code{thresholds}: A numeric vector specifying the detection thresholds for plotting. Default is \code{10^seq(log10(3e-2), log10(1), length = 10)}.
#'     \item \code{min_prevalence}: A numeric value specifying the minimum prevalence threshold for core microbiome. Default is 0.2.
#'     \item \code{taxa_order}: A character string indicating whether to order taxa by "ascending" or "descending" abundance. Default is "descending".
#'   }
#' @param output_core_csv A character string specifying the path to save the core microbiome subset as a CSV file. Default is \code{NULL}, meaning no CSV file is saved.
#' @param output_core_rds A character string specifying the path to save the core microbiome subset as an RDS file. Default is \code{NULL}, meaning no RDS file is saved.
#' @return A \code{ggplot2} object representing the core microbiome prevalence heatmap.
#' @importFrom phyloseq tax_glom prune_taxa subset_taxa transform_sample_counts taxa_names psmelt
#' @importFrom microbiomeutilities format_to_besthit
#' @importFrom microbiome plot_core
#' @importFrom ggplot2 ggplot theme xlab element_text element_blank element_line element_rect scale_x_discrete
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils write.csv
#' @examples
#' \dontrun{
#' # Example usage:
#' custom_detections <- list(
#'   prevalences = seq(0.05, 1, 0.01),  # Custom prevalences
#'   thresholds = 10^seq(log10(0.05), log10(1), length = 8),  # detection thresholds
#'   min_prevalence = 0.4,  # Custom minimum prevalence
#'   taxa_order = "ascending"  # Order taxa by ascending/desc abundance
#' )
#'
#' # Create and save the plot
#' plot_result <- plot_core_microbiome_custom(physeq = ps, 
#'                                            detections = custom_detections, 
#'                                            taxrank = "Family", 
#'                                            output_core_rds = "core_microbiome.rds", 
#'                                            output_core_csv = "core_microbiome.csv")
#'
#' # Print the plot
#' print(plot_result)
#' }
#' @export
plot_core_microbiome_custom <- function(physeq = NULL, 
                                        taxrank = "Genus", 
                                        select_taxa = NULL, 
                                        detections = list(prevalences = seq(0.03, 1, 0.01), 
                                                          thresholds = 10^seq(log10(3e-2), log10(1), length = 10),
                                                          min_prevalence = 0.2,  # Default 20%
                                                          taxa_order = "descending"),  # Default order
                                        output_core_csv = NULL, 
                                        output_core_rds = NULL) {
  
  # Check for required input
  if (is.null(physeq)) {
    stop("Error: 'physeq' argument is required.")
  }
  
  # Provide default values for custom_detections if they are not specified by the user
  detections <- utils::modifyList(list(
    prevalences = seq(0.03, 1, 0.01),  # Default prevalence thresholds
    thresholds = 10^seq(log10(3e-2), log10(1), length = 10),  # Default detection thresholds
    min_prevalence = 0.2,  # Default minimum prevalence (20%)
    taxa_order = "descending"  # Default order of taxa
  ), detections)
  
  # Extract components from the detections list
  prevalences <- detections$prevalences
  thresholds <- detections$thresholds
  min_prevalence <- detections$min_prevalence
  taxa_order <- detections$taxa_order  # Order taxa by ascending/descending
  
  # Glom taxa at specified taxonomic rank
  glom_phy <- phyloseq::tax_glom(physeq, taxrank = taxrank)
  
  # Prune taxa if specific taxa are selected
  if (!is.null(select_taxa)) {
    glom_phy <- phyloseq::prune_taxa(select_taxa, glom_phy)
  }
  
  # Exclude the species level/Sp
  glom_phy <- phyloseq::subset_taxa(glom_phy, select = -Species)
  
  # Transform counts to relative abundance/"plot_core of microbiome requirement"
  glom_phy <- phyloseq::transform_sample_counts(glom_phy, function(x) 100 * x / sum(x))
  
  # Rename taxa with ASV prefix
  phyloseq::taxa_names(glom_phy) <- paste0("ASV", seq(phyloseq::ntaxa(glom_phy)))
  
  # Format phyloseq obj for besthit, its more informative
  phy_rel_f <- microbiomeutilities::format_to_besthit(glom_phy)
  
  # Filter out rows with NA values in critical columns
  pm_core <- phyloseq::psmelt(phy_rel_f) %>%
    dplyr::filter(!is.na(Abundance), !is.na(OTU), !is.na(Sample))
  
  # Filter based on detection thresholds and min_prevalence
  pm_core <- dplyr::filter(pm_core, Abundance >= min(thresholds) & Abundance <= max(thresholds)) %>%  # Filter taxa based on detection thresholds
    dplyr::group_by(OTU) %>%
    dplyr::filter(mean(Abundance > 0) >= min_prevalence)  # Filter taxa based on min_prevalence
  
  # Order taxa by abundance or other criteria/dplyr
  if (taxa_order == "descending") {
    pm_core <- dplyr::arrange(pm_core, dplyr::desc(Abundance))  # Order by descending abundance
  } else {
    pm_core <- dplyr::arrange(pm_core, Abundance)  # Order by ascending abundance
  }
  
  # Check if there's any taxa left after processing
  if (nrow(pm_core) == 0) {
    stop("Error: No taxa remaining after filtering based on detection thresholds and prevalence.")
  }
  
  # Save core microbiome as CSV file if output_core_csv is provided
  if (!is.null(output_core_csv)) {
    utils::write.csv(pm_core, output_core_csv, row.names = FALSE)
    cat("Core microbiome saved as CSV to:", output_core_csv, "\n")
  }
  
  # Save core microbiome as RDS file if output_core_rds is provided
  if (!is.null(output_core_rds)) {
    saveRDS(phy_rel_f, file = output_core_rds)
    cat("Core microbiome saved as RDS to:", output_core_rds, "\n")
  }
  
  # Suppress messages and warnings from plot_core
  p.core <- tryCatch({
    suppressWarnings(suppressMessages(
      microbiome::plot_core(phy_rel_f, 
                            plot.type = "heatmap", 
                            colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
                            prevalences = prevalences, 
                            detections = thresholds,  # Use the custom detection thresholds
                            min.prevalence = min_prevalence)
    )) +  # Use min.prevalence to zoom in
      ggplot2::xlab("Detection Threshold (Relative Abundance)") +  # X-axis label
      ggplot2::scale_x_discrete(labels = function(x) sprintf("%.2f", as.numeric(x))) +  # Format x-axis labels to 2 decimal places
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),  # Remove background
        panel.border = ggplot2::element_blank(),  # Remove border
        panel.grid = ggplot2::element_blank(),  # Remove grid lines
        axis.line = ggplot2::element_line(colour = "black"),  # Keep axis lines
        axis.title.x = ggplot2::element_text(size = 12, face = "bold"),  # Bold X-axis title
        axis.title.y = ggplot2::element_blank(),  # Remove Y-axis title
        axis.text.x = ggplot2::element_text(size = 11, face = "bold", angle = 25, hjust = 1),  # Bold X-axis text rotated 45 degrees
        axis.text.y = ggplot2::element_text(size = 11, face = "bold"),  # Bold Y-axis text
        legend.text = ggplot2::element_text(size = 12, face = "bold"),  # Bold legend text
        legend.title = ggplot2::element_text(size = 13, face = "bold"),  # Bold legend title
        legend.key.size = ggplot2::unit(1, 'cm')  # Adjust legend key size
      )
  }, error = function(e) {
    message("Plotting failed: ", e)
    NULL
  })
  
  return(p.core)
}


#custom_detections <- list(
#  prevalences = seq(0.03, 1, 0.01),  # Custom prevalences
#  thresholds = 10^seq(log10(0.03), log10(1), length = 10),  # Custom thresholds
#  min_prevalence = 0.3,  
#  taxa_order = "ascending") # Order taxa by ascending/descending abundance

# Create and save the plot
#plot_result <- plot_core_microbiome_custom(
#  physeq = ps, 
#  detections = custom_detections, 
# taxrank = "Family", 
#output_core_rds = "core_microbiome.rds", 
#output_core_csv = "core_microbiome.csv")

# Print the plot
#print(plot_result)
