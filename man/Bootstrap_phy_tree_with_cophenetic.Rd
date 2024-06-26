% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap_phy_tree_with_cophenetic.R
\name{Bootstrap_phy_tree_with_cophenetic}
\alias{Bootstrap_phy_tree_with_cophenetic}
\title{Create and Plot a Phylogenetic Tree with Bootstrap Values and Cophenetic Distances}
\usage{
Bootstrap_phy_tree_with_cophenetic(
  physeq_object,
  output_file = "tree_with_bootstrap_and_cophenetic.png",
  bootstrap_replicates = 100
)
}
\arguments{
\item{physeq_object}{A phyloseq object containing the phylogenetic data.}

\item{output_file}{A character string specifying the output file path for the plot. Default is "tree_with_bootstrap_and_cophenetic.png".}

\item{bootstrap_replicates}{An integer specifying the number of bootstrap replicates. Default is 100.}
}
\value{
NULL. The function saves the plot to the specified output file.
}
\description{
This function generates a phylogenetic tree with bootstrap values and cophenetic distances.
It performs multiple sequence alignment, computes bootstrap values, and plots the tree.
}
\examples{
# Bootstrap and plot a phylogenetic tree
Bootstrap_phy_tree_with_cophenetic(physeq_object = Tetragenococcus, 
output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 500)
}
