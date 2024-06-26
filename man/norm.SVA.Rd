% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_set.R
\name{norm.SVA}
\alias{norm.SVA}
\title{Surrogate Variable Analysis for Sequencing Data (SVA)}
\usage{
norm.SVA(raw, groups)
}
\arguments{
\item{raw}{A matrix of raw counts.}

\item{groups}{A vector of group labels.}
}
\value{
A list containing the normalized data and adjustment factor.
}
\description{
Identifies and removes unwanted variation in sequencing data.
}
\references{
SVA package.
}
