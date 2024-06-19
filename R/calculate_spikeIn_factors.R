#' Calculate Scaling Factors for Spiked Species in Phyloseq Object
#'
#' This function calculates the scaling factors for specified spiked species in a phyloseq object.
#' It excludes the spiked species, merges the spiked species into one ASV, calculates the scaling factors,
#' and saves the results as RDS, CSV, and DOCX files.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param spiked_cells A numeric value specifying the number of spiked cells.
#' @param merged_spiked_species A character vector of spiked species to merge in the phyloseq object.
#' @param output_path A character string specifying the path to save the DOCX output file. Default is NULL, which saves the report as "spikeIn_factors_summary.docx".
#' @return A list containing the scaling factors and other related outputs.
#' @examples
#' # Example usage:
#' # Define the spiked species
#' merged_spiked_species <- c("Tetragenococcus_halophilus")
#'
#' # Calculate scaling factors and generate the report
#' result <- calculate_spikeIn_factors(spiked_16S, 1874, merged_spiked_species)
#'
#' # Access the results
#' scaling_factors <- result$scaling_factors
#' physeq_no_spiked <- result$physeq_no_spiked
#' spiked_16S_total_reads <- result$spiked_16S_total_reads
#' spiked_species <- result$spiked_species
#' spiked_species_merged <- result$spiked_species_merged
#' spiked_species_reads <- result$spiked_species_reads
#' @export
calculate_spikeIn_factors <- function(physeq, spiked_cells, merged_spiked_species, output_path = NULL) {
  output_prefix <- "spikeIn_factors"  # Define the fixed output directory to have the files in one folder
  
  cat("Identifying spiked species...\n")
  
  # Check if the output_prefix directory exists, create it if necessary
  if (!dir.exists(output_prefix)) {
    if (!dir.create(output_prefix)) {
      stop("Failed to create directory:", output_prefix)
    }
  }
  
  # Exclude spiked species and save it
  physeq_no_spiked <- subset_taxa(physeq, !Species %in% merged_spiked_species)
  saveRDS(physeq_no_spiked, file.path(output_prefix, "Physeq_NO_Spiked_Species.rds"))
  cat("Saved file:", file.path(output_prefix, "Physeq_NO_Spiked_Species.rds\n"))
  
  # Calculate total reads
  spiked_16S_total_reads <- data.frame(Sample = rownames(sample_data(physeq)), Total_Reads = sample_sums(otu_table(physeq)))
  write.csv(spiked_16S_total_reads, file.path(output_prefix, "Spiked_16S_Total_Reads.csv"))
  cat("Saved file:", file.path(output_prefix, "Spiked_16S_Total_Reads.csv\n"))
  
  # Subset taxa with specified spiked species
  spiked_species <- subset_taxa(physeq, Species %in% merged_spiked_species)
  saveRDS(spiked_species, file.path(output_prefix, "Spiked_Species.rds"))
  
  # Merge spiked species into one
  spiked_species_merged <- merge_taxa(spiked_species, taxa_names(spiked_species))
  saveRDS(spiked_species_merged, file.path(output_prefix, "Spiked_Species_MergedToOne.rds"))
  cat("Saved file:", file.path(output_prefix, "Spiked_Species_MergedToOne.rds\n"))
  
  # Calculate total reads for each spiked species
  spiked_species_reads <- data.frame(Sample = rownames(sample_data(spiked_species_merged)), Total_Reads = sample_sums(otu_table(spiked_species_merged)))
  write.csv(spiked_species_reads, file.path(output_prefix, "Spiked_Species_Reads.csv"))
  cat("Saved file:", file.path(output_prefix, "Spiked_Species_Reads.csv\n"))
  
  # Check if any negative values are present in the spiked species reads
  if (any(spiked_species_reads$Total_Reads <= 0)) {
    cat("Negative values detected in spiked species reads. Converting to zero.\n")
    spiked_species_reads$Total_Reads[spiked_species_reads$Total_Reads <= 0] <- 0
  }
  
  # Calculate scaling factors for each spiked species
  metadata <- sample_data(physeq)  # Access sample data
  scaling_factors <- ifelse(spiked_species_reads$Total_Reads != 0, {
    ifelse(metadata$spiked_volume == 0.5, {
      spiked_cells / 4 / spiked_species_reads$Total_Reads
    }, ifelse(metadata$spiked_volume == 1, {
      spiked_cells / 2 / spiked_species_reads$Total_Reads
    }, ifelse(metadata$spiked_volume == 2, {
      spiked_cells / spiked_species_reads$Total_Reads
    }, ifelse(metadata$spiked_volume == 3, {
      spiked_cells * 1.5 / spiked_species_reads$Total_Reads
    }, ifelse(metadata$spiked_volume == 4, {
      spiked_cells * 2 / spiked_species_reads$Total_Reads
    }, 0)))))
  }, 0)
  
  # Change zeros in scaling factors to 1
  scaling_factors[scaling_factors == 0] <- 1
  
  # Save the scaling factors as 'spikeIn_factors.csv'
  write.csv(scaling_factors, file.path(output_prefix, "spikeIn_factors.csv"), row.names = FALSE)
  cat("Scaling factors saved as spikeIn_factors.csv.\n")
  
  cat("Process complete. Scaling factors for spiked species saved in the output directory.\n")
  
  # Create a data frame for flextable
  scaling_factor_df <- data.frame(Sample = rownames(spiked_species_reads), Scaling_Factor = scaling_factors)
  
  # Create flextable
  ft <- flextable::flextable(scaling_factor_df)
  ft <- flextable::fontsize(ft, size = 10)
  ft <- flextable::font(ft, part = "all", fontname = "Inconsolata")
  ft <- flextable::color(ft, part = "header", color = "red4")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::italic(ft)
  
  # Set default output directory if none provided
  if (is.null(output_path)) {
    output_path <- "spikeIn_factors_summary.docx"
  }
  
  # Save the flextable as a Word document
  flextable::save_as_docx(ft, path = output_path)
  
  # Print a message 
  cat("Table saved in docx format:", output_path, "\n")
  
  # Return the scaling factors along with other outputs
  return(list(scaling_factors = scaling_factors,
              physeq_no_spiked = physeq_no_spiked,
              spiked_16S_total_reads = spiked_16S_total_reads,
              spiked_species = spiked_species,
              spiked_species_merged = spiked_species_merged,
              spiked_species_reads = spiked_species_reads))
}

# Example usage:
#merged_spiked_species <- c("Tetragenococcus_halophilus")
#result <- calculate_spikeIn_factors(physeqASV16, 1874, merged_spiked_species)


#merged_spiked_species <- "Dekkera_bruxellensis"
#spiked_cells<- 733
#result <- calculate_spikeIn_factors(physeq_ITS, 733, merged_spiked_species)

# Access the results
# scaling_factors <- result$scaling_factors
# physeq_no_spiked <- result$physeq_no_spiked
# spiked_16S_total_reads <- result$spiked_16S_total_reads
# spiked_species <- result$spiked_species
# spiked_species_merged <- result$spiked_species_merged
# spiked_species_reads <- result$spiked_species_reads
