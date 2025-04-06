#' @title Plot Core Microbiome Prevalence Heatmap (Phyloseq & TSE Compatible)
#'
#' @description This function generates a prevalence heatmap of the core microbiome at a specified taxonomic rank.
#' It allows users to pass custom detection thresholds, prevalence thresholds, and a minimum prevalence
#' filter, and it provides an option to order taxa either in ascending or descending abundance.
#' The plot displays a heatmap showing detection thresholds at different prevalence levels.
#'
#' @param obj A \code{phyloseq} or \code{TreeSummarizedExperiment} (TSE) object containing microbiome data.
#' @param taxrank A character string specifying the taxonomic rank to glom taxa. Default is "Genus".
#' @param select_taxa A character vector of taxa to select. Default is \code{NULL}, meaning no specific taxa are selected.
#' @param detections A list with the following elements:
#'   \itemize{
#'     \item \code{prevalences}: A numeric vector specifying the prevalence thresholds for plotting. Default is \code{seq(0.03, 1, 0.01)}.
#'     \item \code{thresholds}: A numeric vector specifying the detection thresholds for plotting. Default is \code{10^seq(log10(3e-2), log10(1), length = 10)}.
#'     \item \code{min_prevalence}: A numeric value specifying the minimum prevalence threshold for core microbiome. Default is 0.2.
#'     \item \code{taxa_order}: A character string indicating whether to order taxa by "ascending" or "descending" abundance. Default is "descending".
#'   }
#' @param output_core_csv Path to save the core microbiome subset as a CSV file.
#'   Default is \code{NULL}. To avoid writing to the working directory during examples,
#'   use \code{file.path(tempdir(), "core_microbiome.csv")}.
#'
#' @param output_core_rds Path to save the core microbiome subset as an RDS file.
#'   Default is \code{NULL}. Use \code{file.path(tempdir(), "core_microbiome.rds")}
#'   to write to a temporary directory during examples or checks.
#' @return A \code{ggplot2} object representing the core microbiome prevalence heatmap.
#' @source Uses microbiome::plot_core() for core heatmap visualization
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Define custom detection parameters
#'   custom_detections <- list(
#'     prevalences = seq(0.03, 1, 0.01),
#'     thresholds = 10^seq(log10(0.03), log10(1), length = 10),
#'     min_prevalence = 0.3,
#'     taxa_order = "ascending"
#'   )
#'
#'   # Define output paths using tempdir()
#'   core_csv <- file.path(tempdir(), "core_microbiome.csv")
#'   core_rds <- file.path(tempdir(), "core_microbiome.rds")
#'
#'   # Generate a core microbiome plot with custom settings
#'   plot_result <- plot_core_microbiome_custom(
#'     obj = physeq_16SOTU,
#'     detections = custom_detections,
#'     taxrank = "Genus",
#'     output_core_rds = core_rds,
#'     output_core_csv = core_csv
#'   )
#'
#'   # Print the resulting plot
#'   print(plot_result)
#' }
#' @importFrom phyloseq tax_glom prune_taxa subset_taxa transform_sample_counts taxa_names psmelt
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom microbiome plot_core
#' @importFrom ggplot2 ggplot theme xlab element_text element_blank element_line element_rect scale_x_discrete
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils write.csv
#'
#' @export
plot_core_microbiome_custom <- function(obj,
                                        taxrank = "Genus",
                                        select_taxa = NULL,
                                        detections = list(
                                          prevalences = seq(0.03, 1, 0.01),
                                          thresholds = 10^seq(log10(3e-2), log10(1), length = 10),
                                          min_prevalence = 0.2,
                                          taxa_order = "descending"
                                        ),
                                        output_core_csv = NULL,
                                        output_core_rds = NULL) {
  # Validate input
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Error: 'obj' must be a phyloseq or TreeSummarizedExperiment (TSE) object.")
  }

  # Convert TSE to phyloseq using the internal function from your package
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj) # Using your package's function
  }

  # Check if tax_table exists before proceeding
  if (is.null(tax_table(obj, errorIfNULL = FALSE))) {
    stop("Error: The phyloseq object must contain a taxonomy table for tax_glom().")
  }

  # Aggregate taxa at the selected rank
  glommed_obj <- tax_glom(obj, taxrank = taxrank)

  # Prune taxa if select_taxa is provided
  if (!is.null(select_taxa)) {
    glommed_obj <- prune_taxa(select_taxa, glommed_obj)
  }

  # Replace OTU hashes with names from the selected taxonomic rank
  if (!is.null(tax_table(glommed_obj, errorIfNULL = FALSE))) {
    tax <- as.data.frame(tax_table(glommed_obj))
    rank_column <- tax[[taxrank]]
    names(rank_column) <- rownames(tax)
    taxa_names(glommed_obj) <- make.unique(as.character(rank_column))
  }

  # Transform counts to relative abundance after renaming taxa
  rel_abund_obj <- transform_sample_counts(glommed_obj, function(x) 100 * x / sum(x))

  # Convert to long format
  pm_core <- psmelt(rel_abund_obj) |>
    dplyr::filter(!is.na(Abundance), !is.na(OTU), !is.na(Sample)) |>
    dplyr::group_by(OTU) |>
    dplyr::filter(mean(Abundance > 0) >= detections$min_prevalence)

  # Order taxa by abundance if requested
  pm_core <- if (detections$taxa_order == "descending") {
    dplyr::arrange(pm_core, dplyr::desc(Abundance))
  } else {
    dplyr::arrange(pm_core, Abundance)
  }

  # Ensure at least one taxa is left after filtering
  if (nrow(pm_core) == 0) {
    stop("Error: No taxa remaining after filtering based on detection thresholds and prevalence.")
  }

  # Save results if output paths are provided
  if (!is.null(output_core_csv)) {
    utils::write.csv(pm_core, output_core_csv, row.names = FALSE)
    message("Core microbiome saved as CSV to: ", output_core_csv)
  }
  if (!is.null(output_core_rds)) {
    saveRDS(rel_abund_obj, file = output_core_rds)
    message("Core microbiome saved as RDS to: ", output_core_rds)
  }

  p.core <- withCallingHandlers(
    microbiome::plot_core(
      rel_abund_obj,
      plot.type = "heatmap",
      colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
      prevalences = detections$prevalences,
      detections = detections$thresholds,
      min.prevalence = detections$min_prevalence
    ),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning")
  ) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Detection Threshold (Relative Abundance)") +
    ggplot2::ylab(NULL) +
    ggplot2::scale_x_discrete(labels = function(x) sprintf("%.2f", as.numeric(x))) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11, angle = 25, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.key.size = ggplot2::unit(1, "cm"),
      panel.border = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "black", size = 0.6),
      axis.line.y = ggplot2::element_line(colour = "black", size = 0.6)
    )
  return(p.core)
}


#  custom_detections <- list(
#  prevalences = seq(0.03, 1, 0.01),  # Custom prevalences
#  thresholds = 10^seq(log10(0.03), log10(1), length = 10),  # Custom thresholds
#  min_prevalence = 0.3,
# taxa_order = "ascending") # Order taxa by ascending/descending abundance

# Create and save the plot
# plot_result <- plot_core_microbiome_custom(
# obj = phy_M_M_ps,
# detections = custom_detections,
# taxrank = "Genus",
# output_core_rds = "core_microbiome.rds",
# output_core_csv = "core_microbiome.csv")

# Print the plot
# print(plot_result)
