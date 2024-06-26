% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_set.R
\name{norm.clr}
\alias{norm.clr}
\title{CLR Normalization (Centered Log-Ratio Transformation)}
\usage{
norm.clr(ps, feature_category = c("all", "iqlr", "zero", "lvha"))
}
\arguments{
\item{ps}{A phyloseq object.}

\item{feature_category}{A character string specifying the feature category ("all", "iqlr", "zero", "lvha").}
}
\value{
A list containing the normalized data and scaling factor.
}
\description{
Computes log-ratios relative to the geometric mean of all features.
}
\references{
Aitchison, 1986.
}
