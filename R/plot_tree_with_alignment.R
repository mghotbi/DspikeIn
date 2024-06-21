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
#' # Plot the phylogenetic tree with multiple sequence alignment
#' plot_tree_with_alignment(Tetragenococcus, output_prefix = "tree_alignment", width = 15, height = 15)
#' @export
plot_tree_with_alignment <- function(physeq, output_prefix = "physeq_tree_alignment", width = 15, height = 15) {
  # Check if physeq is a valid phyloseq object
  if (!is(physeq, "phyloseq")) {
    stop("Input 'physeq' must be a valid phyloseq object.")
  }
  
  # Ensure the phy_tree is correctly extracted
  tree <- phy_tree(physeq)
  
  # if tree has branch lengths
  if (is.null(tree$edge.length) || length(tree$edge.length) == 0) {
    stop("Tree has no branch lengths")
  }
  
  # Extract reference sequences from the phyloseq obj
  ref_sequences <- refseq(physeq)
  
  # Perform multiple sequence alignment
  alignment <- DECIPHER::AlignSeqs(ref_sequences, anchor = NA)
  
  # Write aligned sequences to a temporary FASTA file
  temp_fasta <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(alignment, temp_fasta)
  
  # Generate the ggtree plot
  tree_plot <- ggtree(tree, layout = "fan", branch.length = "branch.length", ladderize = TRUE, 
                      right = TRUE, root.position = 0) + theme_tree2()
  
  # Add the multiple sequence alignment plot
  ms <- msaplot(p = tree_plot, fasta = temp_fasta)
  
  # Customize the plot appearance
  final_plot <- ms + scale_color_continuous(low = 'gray88', high = 'red4') + 
    geom_tiplab(aes(color = branch.length), align = TRUE, size = 3) + 
    theme(legend.position = "right", legend.key.size = unit(3, "lines")) + 
    geom_tippoint(color = "#BB650B", shape = "*", size = 3) + theme_tree2()
  
  print(final_plot)
  
  # Save the plots
  png_filename <- paste0(output_prefix, ".png")
  pdf_filename <- paste0(output_prefix, ".pdf")
  
  ggsave(png_filename, plot = final_plot, width = width, height = height, units = "in")
  ggsave(pdf_filename, plot = final_plot, width = width, height = height, units = "in")
  
  cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
}

# Plot the phylogenetic tree with multiple sequence alignment
#plot_tree_with_alignment(Tetragenococcus, output_prefix = "tree_alignment", width = 15, height = 15)
