% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization_set.R
\name{norm.tss}
\alias{norm.tss}
\title{TSS Normalization (Total Sum Scaling)}
\usage{
norm.tss(object)
}
\arguments{
\item{object}{A phyloseq object.}
}
\value{
A phyloseq object with TSS normalized counts.
}
\description{
Converts the feature table into relative abundance by dividing the total reads of each sample.
}
\references{
Bolstad et al., 2003.
}
