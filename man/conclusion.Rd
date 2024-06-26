% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conclusion.R
\name{conclusion}
\alias{conclusion}
\title{Generate Spike Success Report and Summary Statistics}
\usage{
conclusion(
  physeq,
  merged_spiked_species,
  max_passed_range = 11,
  output_path = NULL
)
}
\arguments{
\item{physeq}{A phyloseq object containing the microbial data.}

\item{merged_spiked_species}{A character vector of spiked species to check in the phyloseq object.}

\item{max_passed_range}{A numeric value specifying the maximum acceptable spike percentage. Default is 11.}

\item{output_path}{A character string specifying the path to save the output files. Default is NULL, which saves the report as "spike_success_report.docx".}
}
\value{
A data frame containing the summary statistics of the spike success report.
}
\description{
This function calculates the spike success report for a given phyloseq object, 
summarizes the results, and saves the summary statistics in both DOCX and CSV formats.
It uses the `flextable` package to create and format the tables in the DOCX file.
}
\examples{
# Define the parameters
merged_spiked_species <- c("Dekkera_bruxellensis")
max_passed_range <- 11
output_path <- "spike_success_report.docx"

# Subset the phyloseq object to exclude blanks
physeq_ITS_adj_scaled_perc <- subset_samples(physeq_ITS_OTU, sample_or_blank != "blank")

# Generate the spike success report and summary statistics
conclusion(physeq_ITS_adj_scaled_perc, merged_spiked_species, max_passed_range, output_path)
}
