#' @title Example TreeSummarizedExperiment (TSE) Object with Tree and Reference Sequences
#'
#' @description
#' This dataset contains synthetic microbiome data stored as a \code{TreeSummarizedExperiment} (TSE) object.
#' The TSE structure includes abundance data, sample metadata, taxonomic annotations,
#' a phylogenetic tree linking ASVs, and reference sequences.
#' It is suitable for benchmarking preprocessing workflows involving taxonomic merging or filtering.
#'
#' @format A \code{TreeSummarizedExperiment} object with:
#' \describe{
#'   \item{assays}{Matrix of observed counts or abundances (ASVs × samples).}
#'   \item{rowData}{Taxonomic annotations per ASV.}
#'   \item{colData}{Sample-level metadata.}
#'   \item{rowTree}{A phylogenetic tree object representing ASV relationships.}
#'   \item{referenceSeq}{A \code{DNAStringSet} of DNA sequences matching ASVs.}
#' }
#'
#' @return A \code{TreeSummarizedExperiment} object with full microbiome context (abundance, taxonomy, tree, and reference sequences).
#'
#' @seealso [TreeSummarizedExperiment::TreeSummarizedExperiment()]
#'
#' @usage data(tse)
#'
#' @source Simulated data generated for reproducibility in microbiome pipelines.
#'
#' @examples
#' if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
#'   data(tse, package = "DspikeIn")
#'   tse
#'   SummarizedExperiment::assay(tse)
#'   SummarizedExperiment::colData(tse)
#'   SummarizedExperiment::rowData(tse)
#'   TreeSummarizedExperiment::rowTree(tse)
#'   TreeSummarizedExperiment::referenceSeq(tse)
#' }
"tse"
