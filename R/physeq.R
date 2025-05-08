#' @title General Example Phyloseq Object
#'
#' @description
#' This is a general-purpose \code{phyloseq} object containing microbiome data,
#' simulated for testing and demonstration of microbiome workflows.
#' It contains OTU abundances, taxonomic annotations, and sample metadata.
#'
#' @format A \code{phyloseq} object with:
#' \describe{
#'   \item{otu_table}{Operational Taxonomic Unit (OTU) abundance matrix.}
#'   \item{tax_table}{Taxonomic classification of OTUs.}
#'   \item{sample_data}{Metadata associated with the samples.}
#' }
#'
#' @return A \code{phyloseq} object containing OTU abundances, taxonomy, and sample metadata.
#'
#' @usage data(physeq)
#'
#' @source Simulated dataset for microbiome analysis.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data(physeq)
#'   physeq
#'   phyloseq::sample_names(physeq)
#'   phyloseq::taxa_names(physeq)
#'   summary(physeq)
#' }
"physeq"
