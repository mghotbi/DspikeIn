#' @title Example TreeSummarizedExperiment (TSE) Object
#'
#' @description
#' This dataset contains simulated microbiome data stored as a \code{TreeSummarizedExperiment} (TSE) object.
#' The TSE structure combines abundance data, sample metadata, and optional feature annotations,
#' suitable for downstream statistical analysis of microbiome datasets.
#'
#' @format A \code{TreeSummarizedExperiment} object with:
#' \describe{
#'   \item{assays}{Matrix of observed counts or abundances.}
#'   \item{colData}{Sample metadata table.}
#'   \item{rowData}{Taxonomic annotations or feature metadata (if available).}
#' }
#'
#' @return A \code{TreeSummarizedExperiment} object containing simulated microbiome data.
#'
#' @seealso [TreeSummarizedExperiment::TreeSummarizedExperiment()]
#' @usage data(tse)
#'
#' @source Example dataset for microbiome analysis.
#'
#' @examples
#' data(tse)
#' tse
#' assay(tse)
#' colData(tse)
#' rowData(tse)
"tse"
