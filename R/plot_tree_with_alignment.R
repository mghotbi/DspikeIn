#' Plot a Phylogenetic Tree with Multiple Sequence Alignment
#'
#' This function generates a phylogenetic tree plot with an additional plot of the multiple sequence alignment.
#'
#' @param physeq A phyloseq object containing the phylogenetic data.
#' @param output_prefix A character string specifying the prefix for the output files. Default is "physeq_tree_alignment".
#' @param width A numeric value specifying the width of the output plot. Default is 15.
#' @param height A numeric value specifying the height of the output plot. Default is 15.
#' @return NULL. The function saves the plots to the specified output files.
#' @examples
#' if (interactive()) {
#'   # Plot the phylogenetic tree with multiple sequence alignment
#'   plot_tree_with_alignment(Tetragenococcus, 
#'   output_prefix = "tree_alignment", width = 15, height = 15)
#' }
#' @importFrom phyloseq phy_tree refseq
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings writeXStringSet
#' @importFrom ggtree ggtree theme_tree2 msaplot geom_tiplab geom_tippoint
#' @importFrom ggplot2 scale_color_continuous theme ggsave
#' @importFrom grid unit
#' @export
plot_tree_with_alignment <- function(physeq, output_prefix = "physeq_tree_alignment", width = 15, height = 15) {
  # Check if physeq is a valid phyloseq object using inherits()
  if (!inherits(physeq, "phyloseq")) {
    stop("Input 'physeq' must be a valid phyloseq object.")
  }
  
  # Extract the phylogenetic tree from the phyloseq object
  tree <- phyloseq::phy_tree(physeq)
  
  # Check if the tree has branch lengths
  if (is.null(tree$edge.length) || length(tree$edge.length) == 0) {
    stop("Tree has no branch lengths")
  }
  
  # Extract reference sequences from the phyloseq object
  ref_sequences <- phyloseq::refseq(physeq)
  
  # Perform multiple sequence alignment
  alignment <- DECIPHER::AlignSeqs(ref_sequences, anchor = NA)
  
  # Write aligned sequences to a temporary FASTA file
  temp_fasta <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(alignment, temp_fasta)
  
  # Generate the ggtree plot
  tree_plot <- ggtree::ggtree(tree, layout = "fan", branch.length = "branch.length", ladderize = TRUE, 
                              right = TRUE, root.position = 0) + ggtree::theme_tree2()
  
  # Add the multiple sequence alignment plot
  ms <- ggtree::msaplot(p = tree_plot, fasta = temp_fasta)
  
  # Customize the plot appearance
  final_plot <- ms + ggplot2::scale_color_continuous(low = 'gray88', high = 'red4') + 
    ggtree::geom_tiplab(ggplot2::aes(color = branch.length), align = TRUE, size = 3) + 
    ggplot2::theme(legend.position = "right", legend.key.size = grid::unit(3, "lines")) + 
    ggtree::geom_tippoint(color = "#BB650B", shape = "*", size = 3) + ggtree::theme_tree2()
  
  print(final_plot)
  
  # Save the plots
  png_filename <- paste0(output_prefix, ".png")
  pdf_filename <- paste0(output_prefix, ".pdf")
  
  ggplot2::ggsave(png_filename, plot = final_plot, width = width, height = height, units = "in")
  ggplot2::ggsave(pdf_filename, plot = final_plot, width = width, height = height, units = "in")
  
  cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
}

# Example usage:
# plot_tree_with_alignment(Tetragenococcus, 
# output_prefix = "tree_alignment", width = 15, height = 15)


