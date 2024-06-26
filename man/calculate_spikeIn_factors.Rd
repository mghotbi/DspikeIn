% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_spikeIn_factors.R
\name{calculate_spikeIn_factors}
\alias{calculate_spikeIn_factors}
\title{Calculate Scaling Factors for Spiked Species in Phyloseq Object}
\usage{
calculate_spikeIn_factors(
  physeq,
  spiked_cells,
  merged_spiked_species,
  output_path = NULL
)
}
\arguments{
\item{physeq}{A phyloseq object containing the microbial data.}

\item{spiked_cells}{A numeric value specifying the number of spiked cells.}

\item{merged_spiked_species}{A character vector of spiked species to merge in the phyloseq object.}

\item{output_path}{A character string specifying the path to save the DOCX output file. Default is NULL, which saves the report as "spikeIn_factors_summary.docx".}
}
\value{
A list containing the scaling factors and other related outputs.
}
\description{
This function calculates the scaling factors for specified spiked species in a phyloseq object.
It excludes the spiked species, merges the spiked species into one ASV, calculates the scaling factors,
and saves the results as RDS, CSV, and DOCX files.
}
\details{
Note: The metadata in the phyloseq object must contain a column named 'spiked.volume' with the appropriate values.
}
\examples{
# Example usage:
# Define the spiked species
merged_spiked_species <- c("Tetragenococcus_halophilus")

# Calculate scaling factors and generate the report
result <- calculate_spikeIn_factors(spiked_16S, 1874, merged_spiked_species)

# Access the results
scaling_factors <- result$scaling_factors
physeq_no_spiked <- result$physeq_no_spiked
spiked_16S_total_reads <- result$spiked_16S_total_reads
spiked_species <- result$spiked_species
spiked_species_merged <- result$spiked_species_merged
spiked_species_reads <- result$spiked_species_reads
}
