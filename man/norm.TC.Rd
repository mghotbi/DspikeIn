% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_set.R
\name{norm.TC}
\alias{norm.TC}
\title{TC Normalization (Total Count Scaling)}
\usage{
norm.TC(raw, groups)
}
\arguments{
\item{raw}{A matrix of raw counts.}

\item{groups}{A vector of group labels.}
}
\value{
A list containing the normalized data and scaling factor.
}
\description{
A simple scaling normalization method dividing by the library size.
}
