#' @title Example Phyloseq Object with Tree and Reference Sequences
#'
#' @description
#' A general-purpose \code{phyloseq} object for microbiome method development and testing.
#' This synthetic dataset includes OTU abundances, taxonomic classifications, sample metadata,
#' a phylogenetic tree, and reference sequences.
#'
#' @format A \code{phyloseq} object with:
#' \describe{
<<<<<<< HEAD
#'   \item{otu_table}{OTU abundance matrix (ASVs × samples).}
=======
#'   \item{otu_table}{OTU abundance matrix (ASVs  samples).}
>>>>>>> bioc-rev1
#'   \item{tax_table}{Taxonomic classification of ASVs.}
#'   \item{sample_data}{Sample-level metadata.}
#'   \item{phy_tree}{A rooted phylogenetic tree of ASVs.}
#'   \item{refseq}{DNA reference sequences corresponding to each ASV.}
#' }
#'
#' @return A \code{phyloseq} object with full microbiome data components.
#'
#' @usage data(physeq)
#'
#' @source Synthetic data generated for benchmarking and demonstration.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data(physeq, package = "DspikeIn")
#'   physeq
#'   phyloseq::sample_names(physeq)
#'   phyloseq::taxa_names(physeq)
#'   phyloseq::phy_tree(physeq)
#'   phyloseq::refseq(physeq)
#' }
"physeq"
