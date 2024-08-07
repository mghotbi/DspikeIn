% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microbiome_differential_abundance.R
\name{perform_edgeR}
\alias{perform_edgeR}
\title{Perform Differential Abundance Analysis with edgeR}
\usage{
perform_edgeR(
  ps,
  group_var,
  threshold_percentage,
  threshold_mean_abundance,
  threshold_count,
  threshold_relative_abundance,
  significance_level
)
}
\arguments{
\item{ps}{A phyloseq object.}

\item{group_var}{A string specifying the grouping variable in sample data.}

\item{threshold_percentage}{A numeric value for filtering threshold percentage.}

\item{threshold_mean_abundance}{A numeric value for filtering mean abundance.}

\item{threshold_count}{A numeric value for filtering count.}

\item{threshold_relative_abundance}{A numeric value for filtering relative abundance.}

\item{significance_level}{A numeric value specifying the significance level for filtering significant OTUs.}
}
\value{
A data frame containing the differential abundance results.
}
\description{
Perform Differential Abundance Analysis with edgeR
}
