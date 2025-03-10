#' @title Create and Plot a Phylogenetic Tree with Bootstrap Values and Cophenetic Distances
#' @description Generates a phylogenetic tree with bootstrap values and cophenetic distances.
#' Supports `phyloseq` and `TreeSummarizedExperiment` (TSE) by converting TSE to `phyloseq`.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param output_file A character string specifying the output file path for the plot.
#' @param bootstrap_replicates An integer specifying the number of bootstrap replicates. Default is 100.
#' @return NULL. Saves the plot to the specified output file.
#' @importFrom phyloseq phy_tree refseq
#' @importFrom DECIPHER AlignSeqs
#' @importFrom ape cophenetic.phylo boot.phylo nj nodelabels edgelabels
#' @importFrom phangorn phyDat dist.ml
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Bootstrap and plot a phylogenetic tree
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Subset phyloseq object to include only Tetragenococcus genus
#'   Tetragenococcus <- phyloseq::subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
#'
#'   # Generate phylogenetic tree with bootstrap and cophenetic distance
#'   Bootstrap_phy_tree_with_cophenetic(
#'     obj = Tetragenococcus,
#'     output_file = "tree_with_bootstrap_and_cophenetic.png",
#'     bootstrap_replicates = 500
#'   )
#' }
#' }
#' @export
Bootstrap_phy_tree_with_cophenetic <- function(obj, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 100) {

  # Convert TSE to phyloseq if necessary
  if (inherits(obj, "TreeSummarizedExperiment")) {
    if (!exists("convert_tse_to_phyloseq")) {
      stop("\U0000274C Function 'convert_tse_to_phyloseq()' not found. Please ensure it is available in your package.")
    }
    obj <- convert_tse_to_phyloseq(obj)
    message("Converted TSE object to phyloseq format.")
  }

  # Extract the phylogenetic tree
  tree <- phyloseq::phy_tree(obj)

  # Validate tree structure
  if (is.null(tree$edge.length) || length(tree$edge.length) == 0) {
    stop("Tree has no branch lengths")
  }

  # Extract reference sequences (only if available)
  if (!is.null(phyloseq::refseq(obj))) {
    ref_sequences <- phyloseq::refseq(obj)
  } else {
    stop("\U0000274C Reference sequences are required for alignment but not found.")
  }

  # Perform multiple sequence alignment
  alignment <- DECIPHER::AlignSeqs(ref_sequences, anchor = NA)

  # Convert alignment to DNAStringSet
  aligned_sequences <- as(alignment, "DNAStringSet")

  # Convert to phangorn phyDat format
  phyDat_alignment <- phangorn::phyDat(as(aligned_sequences, "matrix"), type = "DNA")

  # Convert to matrix for bootstrapping
  alignment_matrix <- as.matrix(aligned_sequences)

  # Bootstrap function
  bootstrap_fun <- function(x) ape::nj(phangorn::dist.ml(phangorn::phyDat(x, type = "DNA")))

  # Generate bootstrap values
  bootstrap_values <- ape::boot.phylo(tree, alignment_matrix, FUN = bootstrap_fun, B = bootstrap_replicates)

  # Normalize bootstrap values to percentages
  bootstrap_percentages <- bootstrap_values / bootstrap_replicates * 100

  # Calculate cophenetic distances
  tree_dist <- ape::cophenetic.phylo(tree)
  print("Cophenetic distance matrix:")
  print(tree_dist)  # Print cophenetic distances in console

  # Plot the phylogenetic tree with bootstrap values
  par(mar = c(5, 4, 4, 10))  # Adjust right margin for labels
  plot(tree, main = "Phylogenetic Tree with Bootstrap Values and Cophenetic Distances", cex = 0.9, tip.color = "dodgerblue4")

  # Add bootstrap values as node labels (red)
  ape::nodelabels(round(bootstrap_percentages, 1), cex = 0.9, frame = "none", adj = c(1.0, -0.5), col = "#6A0572")

  # Add cophenetic distances as edge labels (blue)
  ape::edgelabels(round(tree$edge.length, 2), cex = 0.8, col = "#FF5722", frame = "none")

  # Save the plot
  png(output_file, width = 1200, height = 1200)
  par(mar = c(5, 4, 4, 10))
  plot(tree, main = "Phylogenetic Tree with Bootstrap Values and Cophenetic Distances", cex = 0.9, tip.color = "dodgerblue4")
  ape::nodelabels(round(bootstrap_percentages, 1), cex = 0.9, frame = "none", adj = c(1.0, -0.5), col = "#6A0572")
  ape::edgelabels(round(tree$edge.length, 2), cex = 0.8, col = "#FF5722", frame = "none")
  dev.off()

  cat("\U0001F5C2 Phylogenetic tree with bootstrap values and cophenetic distances saved as:", output_file, "\n")
}

# Example:
# Load the data from DspikeIn package
# Tetragenococcus <-subset_taxa(physeq_16SOTU, Genus== "Tetragenococcus")
# Bootstrap_phy_tree_with_cophenetic(obj = Tetragenococcus,
# output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 100)
