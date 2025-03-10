#' @title Plot and Save a Phylogenetic Tree (phyloseq/TSE)
#'
#' @description Plots a phylogenetic tree from either a `phyloseq` or `TreeSummarizedExperiment` object.
#'              If the tree is missing in a TSE, it converts the object to `phyloseq` before plotting.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing phylogenetic data.
#' @param output_prefix A character string specifying the prefix for the output files. Default is `"phy_tree"`.
#' @param width A numeric value specifying the width of the output plot (in inches). Default is `10`.
#' @param height A numeric value specifying the height of the output plot (in inches). Default is `10`.
#' @param layout A character string specifying the tree layout. Options are `"rectangular"`, `"circular"`, `"fan"`, and `"radial"`. Default is `"circular"`.
#' @param tip_labels A logical indicating whether to display tip labels. Default is `TRUE`.
#' @param metadata_column A string specifying a column from `colData(obj)` (TSE) or `sample_data(obj)` (phyloseq) to annotate the tree nodes. Default is `NULL`.
#' @param print_plot A logical indicating whether to print the plot to the console. Default is `TRUE`.
#' @return NULL. Saves the tree plot as PNG and PDF files.
#'
#' @importFrom phyloseq phy_tree sample_data
#' @importFrom TreeSummarizedExperiment rowTree
#' @importFrom SummarizedExperiment colData
#' @importFrom ggtree ggtree geom_tiplab geom_tippoint geom_text theme_tree2
#' @importFrom ggplot2 ggsave aes theme element_text element_blank theme_minimal
#' @importFrom grid unit
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Plot tree for phyloseq object
#'   plot_tree_custom(physeq_16SOTU, output_prefix = "phylo_plot", layout = "circular")
#'
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Example with a TreeSummarizedExperiment object and metadata annotation
#'   plot_tree_custom(
#'     tse_16SOTU,
#'     output_prefix = "tse_plot",
#'     layout = "fan",
#'     metadata_column = "SampleType"
#'   )
#' }
#' }
#' @export
plot_tree_custom <- function(obj, output_prefix = "phy_tree", width = 10, height = 10, layout = "circular",
                             tip_labels = TRUE, metadata_column = NULL, print_plot = TRUE) {
  suppressMessages({

    #  Validate input type
    if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
      stop(" Error: Input must be a 'phyloseq' or 'TreeSummarizedExperiment' object.")
    }

    #  Try extracting tree directly (for TSE)
    if (inherits(obj, "TreeSummarizedExperiment")) {
      tree <- TreeSummarizedExperiment::rowTree(obj)

      # If tree is NULL, try converting to phyloseq
      if (is.null(tree)) {
        message(" No tree found in TSE. Attempting to convert TSE to phyloseq...")
        obj <- convert_tse_to_phyloseq(obj)
        tree <- phyloseq::phy_tree(obj)  # Extract tree from converted phyloseq
      }
    } else {
      #  If already a phyloseq object, extract tree normally
      tree <- phyloseq::phy_tree(obj)
    }

    #  Handle missing tree case
    if (is.null(tree)) {
      warning("No phylogenetic tree found. Skipping tree plotting.")
      return(NULL)
    }

    #  Generate the ggtree plot
    tree_plot <- ggtree::ggtree(tree, layout = layout) +
      ggtree::theme_tree2(margin = 10) +  # Increase spacing to avoid cutoff
      ggplot2::theme(
        legend.position = "bottom",  # Move legend below
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 8),
        plot.margin = grid::unit(c(1, 1, 1, 1), "cm")  # Prevent overflow
      )

    #  Add metadata labels (if specified)
    if (!is.null(metadata_column)) {
      metadata_df <- get_sample_data(obj)

      if (!metadata_column %in% colnames(metadata_df)) {
        warning(" Specified metadata column not found in object. Proceeding without annotation.")
      } else {
        metadata_values <- metadata_df[[metadata_column]]

        if (length(metadata_values) == length(tree$tip.label)) {
          tree_plot <- tree_plot +
            ggtree::geom_text(ggplot2::aes(label = metadata_values), size = 2, hjust = 0, color = "#D00000")
        } else {
          warning(" Metadata length does not match the number of tree tips. Skipping metadata annotation.")
        }
      }
    }

    #  Add tip labels if enabled
    if (tip_labels) {
      tree_plot <- tree_plot +
        ggtree::geom_tiplab(ggplot2::aes(color = branch.length), align = TRUE, size = 1.5)  # Smaller labels
    }

    #  Add tip points for better visualization
    tree_plot <- tree_plot +
      ggtree::geom_tippoint(color = "#BB650B", shape = "*", size = 3)

    #  Print plot to console if `print_plot = TRUE`
    if (print_plot) {
      print(tree_plot)
    }

    #  Define output filenames
    png_filename <- paste0(output_prefix, ".png")
    pdf_filename <- paste0(output_prefix, ".pdf")

    #  Save plots
    ggplot2::ggsave(png_filename, plot = tree_plot, width = width, height = height, units = "in")
    ggplot2::ggsave(pdf_filename, plot = tree_plot, width = width, height = height, units = "in")

    message(" Plots saved as: ", png_filename, " and ", pdf_filename)
  })
}

# Example usage:
# Tetragenococcus<-subset_taxa(physeq_16SOTU,Genus=="Tetragenococcus")
# Tetragenococcus_TSE<-convert_phyloseq_to_tse(Tetragenococcus)

# plot_tree_custom(Tetragenococcus, output_prefix = "p0",
# width = 18, height = 18, layout = "circular")

#  Plot phylogenetic tree normally
# plot_tree_custom(Tetragenococcus_TSE,
# output_prefix = "phylo_plot", layout = "circular")

#  Plot a TreeSummarizedExperiment object
# plot_tree_custom(Tetragenococcus_TSE, output_prefix = "tse_plot", layout = "fan")

#  Plot without tip labels (reduces clutter)
# plot_tree_custom(Tetragenococcus_TSE, output_prefix = "phylo_clean", tip_labels = FALSE)
