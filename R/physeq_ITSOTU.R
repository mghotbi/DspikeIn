#' @title Example Phyloseq Object for ITS OTUs
#'
#' @description
#' This dataset contains an example \code{phyloseq} object representing
#' Internal Transcribed Spacer (ITS) amplicon sequencing data. It includes
#' taxonomic annotations, OTU abundance counts, and associated sample metadata,
#' suitable for downstream analysis of fungal communities.
#'
#' @format A \code{phyloseq} object with:
#' \describe{
#'   \item{otu_table}{Operational Taxonomic Unit (OTU) abundance matrix.}
#'   \item{tax_table}{Taxonomic classification of OTUs.}
#'   \item{sample_data}{Metadata associated with the samples.}
#'   \item{phy_tree}{Phylogenetic tree relating OTUs (if available).}
#' }
#'
#' @return A \code{phyloseq} object containing ITS OTU data with taxonomy, sample metadata, and phylogeny.
#'
#' @usage data(physeq_ITSOTU)
#'
#' @source Internal dataset for microbiome analysis.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data(physeq_ITSOTU)
#'   physeq_ITSOTU
#'   summary(physeq_ITSOTU)
#'   phyloseq::sample_data(physeq_ITSOTU)
#'   phyloseq::taxa_names(physeq_ITSOTU)
#' }
"physeq_ITSOTU"
