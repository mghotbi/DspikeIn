#' @title Example Phyloseq Object for 16S OTUs
#'
#' @description
#' This dataset contains an example \code{phyloseq} object representing
#' 16S rRNA gene amplicon sequencing data. It includes taxonomic assignments,
#' abundance counts, sample metadata, and phylogenetic information.
#' This object is intended for demonstration and testing of microbiome workflows.
#'
#' @format A \code{phyloseq} object with:
#' \describe{
#'   \item{otu_table}{Operational Taxonomic Unit (OTU) abundance matrix.}
#'   \item{tax_table}{Taxonomic classification of OTUs.}
#'   \item{sample_data}{Metadata associated with the samples.}
#'   \item{phy_tree}{Phylogenetic tree relating OTUs (if available).}
#' }
#'
#' @return A \code{phyloseq} object containing 16S OTU data with taxonomy, sample metadata, and phylogeny.
#'
#' @usage data(physeq_16SOTU)
#'
#' @source Internal dataset for microbiome analysis.
#'
#' @examples
#' data(physeq_16SOTU)
#' physeq_16SOTU
#' summary(physeq_16SOTU)
#' sample_names(physeq_16SOTU)
#' taxa_names(physeq_16SOTU)
"physeq_16SOTU"
