% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ridge_Plot_it.R
\name{ridge_plot_it}
\alias{ridge_plot_it}
\title{Generate Ridge Plots to Visualize the Distribution of Abundances of Taxa}
\usage{
ridge_plot_it(physeq, taxrank = "Genus", rarefaction_depth = NULL, top_n = 10)
}
\arguments{
\item{physeq}{A phyloseq object containing the taxonomic and abundance data.}

\item{taxrank}{A character string specifying the taxonomic rank to use for glomming and plotting. Default is "Genus".}

\item{rarefaction_depth}{A numeric value specifying the rarefaction depth. If NULL, it is set to 90% of the minimum sample sum. Default is NULL.}

\item{top_n}{An integer specifying the number of top taxa to include in the plot. Default is 10.}
}
\value{
A ggplot object representing the ridge plot of the distribution of relative abundances of taxa.
}
\description{
This function rarefies the data, performs proportion transformation, gloms and prunes taxa, and generates ridge plots to visualize the distribution of relative abundances of taxa.
}
\examples{
# Example usage:
# ridgeP <- ridge_plot_it(spiked_16S, taxrank = "Family", top_n = 10)
# print(ridgeP) + my_custom_theme()
}
