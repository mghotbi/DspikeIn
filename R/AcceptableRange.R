#' @title Acceptable Range Data
#'
#' @description
#' This dataset provides acceptable ranges for microbial abundance estimates,
#' which can be used for validation and quality control purposes in microbiome analyses.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{min_abundance}{The minimum acceptable abundance threshold.}
#'   \item{max_abundance}{The maximum acceptable abundance threshold.}
#'   \item{taxon}{Taxonomic identifier or group name.}
#' }
#'
#' @return A data frame with acceptable abundance thresholds by taxon.
#'
#' @usage data(AcceptableRange)
#'
#' @source Internal package dataset.
#'
#' @examples
#' data(AcceptableRange)
#' head(AcceptableRange)
#' summary(AcceptableRange)
"AcceptableRange"
