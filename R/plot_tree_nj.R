#' Plot a Neighbor-Joining Phylogenetic Tree with Bootstrap Values
#'
#' This function reads DNA sequences from a FASTA file, constructs a Neighbor-Joining phylogenetic tree,
#' generates bootstrap values, and plots the tree.
#'
#' @param fasta_file A character string specifying the path to the FASTA file containing DNA sequences.
#' @param output_file A character string specifying the output file path for the plot. Default is "neighbor_joining_tree.png".
#' @return NULL. The function saves the plot to the specified output file and displays it.
#' @examples
#' \dontrun{
#' # Plot the Neighbor-Joining tree from the "tetra.fasta" file
#' # The output will be saved as "neighbor_joining_tree_with_bootstrap.png"
#' plot_tree_nj("tetra.fasta", output_file = "neighbor_joining_tree_with_bootstrap.png")
#' }
#' @importFrom ape nj plot.phylo boot.phylo nodelabels
#' @importFrom Biostrings readDNAStringSet
#' @importFrom msa msa
#' @importFrom phangorn phyDat dist.ml
#' @export
plot_tree_nj <- function(fasta_file, output_file = "neighbor_joining_tree.png") {
  suppressMessages({
    # Load necessary libraries
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("Package 'Biostrings' is required but not installed.")
    }
    if (!requireNamespace("msa", quietly = TRUE)) {
      stop("Package 'msa' is required but not installed.")
    }
    if (!requireNamespace("phangorn", quietly = TRUE)) {
      stop("Package 'phangorn' is required but not installed.")
    }
    if (!requireNamespace("ape", quietly = TRUE)) {
      stop("Package 'ape' is required but not installed.")
    }
    
    # Read the DNA sequences from a FASTA file
    sequence <- Biostrings::readDNAStringSet(fasta_file)
    
    # Perform multiple sequence alignment
    my_alignment <- msa::msa(sequence)
    
    # Convert alignment to DNAStringSet
    aligned_sequences <- as(my_alignment, "DNAStringSet")
    
    # Convert DNAStringSet to phyDat format
    phyDat_alignment <- phangorn::phyDat(as(aligned_sequences, "matrix"), type = "DNA")
    
    # Compute distance matrix using maximum likelihood
    distance_matrix <- phangorn::dist.ml(phyDat_alignment)
    
    # Construct the phylogenetic tree using the Neighbor-Joining method
    phylo_tree <- ape::nj(distance_matrix)
    
    # Generate bootstrap values
    set.seed(123)  # For reproducibility
    bootstrap_values <- ape::boot.phylo(phylo_tree, as.matrix(aligned_sequences), 
                                        FUN = function(x) ape::nj(phangorn::dist.ml(phangorn::phyDat(x, type = "DNA"))), B = 100)
    
    # Normalize bootstrap values to percentages
    bootstrap_percentages <- bootstrap_values / 100 * 100
    
    # Plot the phylogenetic tree with bootstrap values
    png(output_file, width = 800, height = 800)
    par(mar = c(5, 4, 4, 10))  # Increase the right margin for node labels
    ape::plot.phylo(phylo_tree, main = "Neighbor Joining Tree with Bootstrap Values", cex = 0.6, tip.color = "dodgerblue4")
    ape::nodelabels(round(bootstrap_percentages, 1), cex = 0.8, frame = "none", adj = c(0.8, -0.1), col = "red4")
    dev.off()
    
    # Display the plot
    ape::plot.phylo(phylo_tree, main = "Neighbor Joining Tree with Bootstrap Values", cex = 0.6, tip.color = "dodgerblue4")
    ape::nodelabels(round(bootstrap_percentages, 1), cex = 0.8, frame = "none", adj = c(0.8, -0.1), col = "red4")
    
    cat("Neighbor Joining tree with bootstrap values saved as:", output_file, "\n")
  })
}

# Example usage:
# Plot Neighbor-Joining tree with bootstrap values
# fasta_path <- "C:/Users/mghotbi/Desktop/DspikeIn_R/Data/Test Data/OTUs/16S/physeqOTU16S/tetra.fasta"
# plot_tree_nj(fasta_path, output_file = "neighbor_joining_tree_with_bootstrap.png")
