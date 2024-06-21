#' Create and Plot a Phylogenetic Tree with Bootstrap Values and Cophenetic Distances
#'
#' This function generates a phylogenetic tree with bootstrap values and cophenetic distances.
#' It performs multiple sequence alignment, computes bootstrap values, and plots the tree.
#'
#' @param physeq_object A phyloseq object containing the phylogenetic data.
#' @param output_file A character string specifying the output file path for the plot. Default is "tree_with_bootstrap_and_cophenetic.png".
#' @param bootstrap_replicates An integer specifying the number of bootstrap replicates. Default is 100.
#' @return NULL. The function saves the plot to the specified output file.
#' @examples
#' # Bootstrap and plot a phylogenetic tree
#' Bootstrap_phy_tree_with_cophenetic(physeq_object = Tetragenococcus, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 500)
#' @export
Bootstrap_phy_tree_with_cophenetic <- function(physeq_object, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 100) {
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("Package 'DECIPHER' is required but not installed.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required but not installed.")
  }
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    stop("Package 'phangorn' is required but not installed.")
  }
  
  library(phyloseq)
  library(DECIPHER)
  library(ape)
  library(phangorn)
  
  # Extract the phylogenetic tree from the phyloseq object
  tree <- phy_tree(physeq_object)
  
  # Check if the tree has branch lengths
  if (is.null(tree$edge.length) || length(tree$edge.length) == 0) {
    stop("Tree has no branch lengths")
  }
  
  # Perform multiple sequence alignment (assuming you have the sequences in your physeq object)
  ref_sequences <- refseq(physeq_object)
  alignment <- DECIPHER::AlignSeqs(ref_sequences, anchor = NA)
  
  # Convert alignment to DNAStringSet
  aligned_sequences <- as(alignment, "DNAStringSet")
  
  # Convert DNAStringSet to phyDat format
  phyDat_alignment <- phyDat(as(aligned_sequences, "matrix"), type = "DNA")
  
  # Convert phyDat_alignment to matrix format for bootstrapping
  alignment_matrix <- as.matrix(aligned_sequences)
  
  # Function for bootstrap
  bootstrap_fun <- function(x) nj(dist.ml(phyDat(x, type = "DNA")))
  
  # Generate bootstrap values
  bootstrap_values <- boot.phylo(tree, alignment_matrix, FUN = bootstrap_fun, B = bootstrap_replicates)
  
  # Normalize bootstrap values to percentages
  bootstrap_percentages <- bootstrap_values / bootstrap_replicates * 100
  
  # Calculate cophenetic distances
  tree_dist <- ape::cophenetic.phylo(tree)
  print(tree_dist)
  
  # Plot the phylogenetic tree with bootstrap values
  par(mar = c(5, 4, 4, 10))  # Increase the right margin for node labels
  plot(tree, main = "Phylogenetic Tree with Bootstrap Values and Cophenetic Distances", cex = 0.9, tip.color = "dodgerblue4")
  nodelabels(round(bootstrap_percentages, 1), cex = 0.9, frame = "none", adj = c(1.0, -0.5), col = "red4")
  
  # Save the plot as PNG
  png(output_file, width = 1200, height = 1200)
  par(mar = c(5, 4, 4, 10))  # Increase the right margin for node labels
  plot(tree, main = "Phylogenetic Tree with Bootstrap Values and Cophenetic Distances", cex = 0.9, tip.color = "dodgerblue4")
  nodelabels(round(bootstrap_percentages, 1), cex = 0.9, frame = "none", adj = c(1.0, -0.5), col = "red4")
  dev.off()
  
  cat("Phylogenetic tree with bootstrap values and cophenetic distances saved as:", output_file, "\n")
}

# Example:
# Bootstrap_phy_tree_with_cophenetic(physeq_object = Tetragenococcus, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 1000)

