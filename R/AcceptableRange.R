#' @title Acceptable Range Data
#'
#' @description
#' This dataset provides reference ranges and sample-level metadata for microbial
#' spike-in performance evaluation. It includes taxonomic annotations and summary
#' statistics used in validation and quality control workflows.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Ecoregion_III}{EPA Level III ecoregion classification.}
#'   \item{Genus}{Genus of the taxon.}
#'   \item{Host_genus}{Host genus from which the sample was derived.}
#'   \item{Percentage}{Observed spike-in percentage.}
#'   \item{Phylum}{Phylum classification of the taxon.}
#'   \item{Range}{Acceptable range category.}
#'   \item{Total_Reads_spiked}{Number of reads matching the spike-in species.}
#'   \item{Total_Reads_total}{Total number of reads per sample.}
#'   \item{X}{Row identifier (optional, may be an index or artifact of data processing).}
#'   \item{mean_abundance}{Mean abundance across all samples.}
#' }
#'
#' @return A data frame of spike-in evaluation metrics and taxonomy annotations.
#'
#' @usage data(AcceptableRange)
#'
#' @source Internal package dataset.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("AcceptableRange", package = "DspikeIn")
#'   head(AcceptableRange)
#'   summary(AcceptableRange$Percentage)
#' }
"AcceptableRange"
