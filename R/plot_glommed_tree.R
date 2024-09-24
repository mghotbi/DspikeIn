#' Plot a Phylogenetic Tree with Glommed OTUs
#'
#' This function gloms the OTUs in the phyloseq object at the specified resolution and plots the tree.
#'
#' @param physeq A phyloseq object containing the phylogenetic data.
#' @param resolution A numeric value specifying the resolution for glomming the OTUs.
#' @param output_prefix A character string specifying the prefix for the output files. Default is "glommed_tree".
#' @param width A numeric value specifying the width of the output plot. Default is 18.
#' @param height A numeric value specifying the height of the output plot. Default is 18.
#' @param print_plot A logical value specifying whether to print the plot. Default is TRUE.
#' @return NULL. The function saves the plot to the specified output file.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Plot the tree with glommed OTUs
#'   plot_glommed_tree(Tetragenococcus, resolution = 0.2, 
#'   output_prefix = "top", width = 18, height = 18)
#' }
#' }
#' @importFrom phyloseq phy_tree
#' @importFrom speedyseq tree_glom
#' @importFrom ggtree ggtree geom_tiplab theme_tree2 geom_text2 
#' @importFrom ggplot2 ggsave aes ggtitle
#' @importFrom tibble as_tibble
#' @export
plot_glommed_tree <- function(physeq, resolution, output_prefix = "glommed_tree", width = 20, height = 18, print_plot = TRUE) {
  
  # Perform tree glomming
  physeq_glommed <- speedyseq::tree_glom(physeq, resolution = resolution)
  
  # Check if tree glomming worked correctly
  if (is.null(physeq_glommed)) {
    stop("Tree glomming failed. Please check the input phyloseq object and resolution.")
  }
  
  # Extract tree from phyloseq object
  phy_tree <- phyloseq::phy_tree(physeq_glommed)
  
  # Convert to a format suitable for ggtree
  tree_data <- tibble::as_tibble(phy_tree)
  
  # Determine if nodes are tips
  tree_data$isTip <- !is.na(match(tree_data$label, phy_tree$tip.label))
  
  # Create ggtree object
  tree_plot <- ggtree::ggtree(phy_tree) +
    ggtree::geom_tiplab(ggplot2::aes(label = label), size = 3) +  # Label tips with ASV names
    ggtree::theme_tree2() +
    ggplot2::ggtitle(paste("Glommed Tree at Resolution", resolution))
  
  # Add bootstrap values if available
  if (!is.null(phy_tree$node.label)) {
    bootstrap_values <- phy_tree$node.label
    internal_nodes <- which(!tree_data$isTip)
    tree_plot <- tree_plot + ggtree::geom_text2(ggplot2::aes(subset = !isTip, label = bootstrap_values[match(node, internal_nodes)]), size = 2, vjust = -0.5)
  } else {
    warning("Bootstrap values not found. Skipping bootstrap labels.")
  }
  
  # Print the plot if required
  if (print_plot) {
    print(tree_plot)
  }
  
  # Save the plot
  ggplot2::ggsave(paste0(output_prefix, "_tree.pdf"), plot = tree_plot, width = width, height = height)
  cat("Glommed tree plot saved to:", paste0(output_prefix, "_tree.pdf"), "\n")
  
  return(NULL)
}

# Example usage:
#plot_glommed_tree(Tetragenococcus, resolution = 0.2, output_prefix = "top", width = 18, height = 18, print_plot = TRUE)
