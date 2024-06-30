% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_set.R
\name{norm.QN}
\alias{norm.QN}
\title{Quantile Normalization (QN)}
\usage{
norm.QN(raw, filter = FALSE)
}
\arguments{
\item{raw}{A matrix of raw counts.}

\item{filter}{Logical, whether to filter low counts.}
}
\value{
A list containing the normalized data.
}
\description{
Adjusts the distribution of counts so that different samples have the same distribution.
}
\references{
Bolstad et al., 2003.
}