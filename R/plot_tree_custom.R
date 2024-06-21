#' Plot and Save a Phylogenetic Tree
#'
#' This function plots a phylogenetic tree from a phyloseq object and saves the plot to a specified file.
#'
#' @param physeq A phyloseq object containing the phylogenetic data.
#' @param output_prefix A character string specifying the prefix for the output files. Default is "phy_tree".
#' @param width A numeric value specifying the width of the output plot. Default is 10.
#' @param height A numeric value specifying the height of the output plot. Default is 10.
#' @param layout A character string specifying the layout of the tree. Options are "rectangular", "circular", "fan", and "radial". Default is "rectangular".
#' @param print_plot A logical value specifying whether to print the plot. Default is TRUE.
#' @return NULL. The function saves the plot to the specified output file.
#' @examples
#' # Plot and save the phylogenetic tree
#' plot_tree_custom(Tetragenococcus, output_prefix = "p0", width = 18, height = 18, layout = "circular")
#' @export
plot_tree_custom <- function(physeq, output_prefix = "phy_tree", width = 10, height = 10, layout = "rectangular", print_plot = TRUE) {
  # Check if physeq is a valid phyloseq object
  if (!inherits(physeq, "phyloseq")) {
    stop("Input 'physeq' must be a valid phyloseq object.")
  }
  
  # Ensure the phy_tree is correctly extracted
  tree <- phyloseq::phy_tree(physeq)
  
  # Create the ggtree plot with specified layout
  tree_plot <- ggtree::ggtree(tree, layout = layout) + ggtree::theme_tree2()
  
  # Customize the plot appearance
  final_plot <- tree_plot + 
    ggtree::geom_tiplab(aes(color = branch.length), align = TRUE, size = 3) + 
    theme(legend.position = "right", legend.key.size = unit(3, "lines")) + 
    ggtree::geom_tippoint(color = "#BB650B", shape = "*", size = 3) + 
    ggtree::theme_tree2()
  
  # Print the plot if required
  if (print_plot) {
    print(final_plot)
  }
  
  # Save the plot
  png_filename <- paste0(output_prefix, ".png")
  pdf_filename <- paste0(output_prefix, ".pdf")
  
  ggplot2::ggsave(png_filename, plot = final_plot, width = width, height = height, units = "in")
  ggplot2::ggsave(pdf_filename, plot = final_plot, width = width, height = height, units = "in")
  
  cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
}

# Example usage:
#plot_tree_custom(Tetragenococcus, output_prefix = "p0", width = 18, height = 18, layout = "circular")
  