pkgname <- "DspikeIn"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DspikeIn')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AcceptableRange")
### * AcceptableRange

flush(stderr()); flush(stdout())

### Name: AcceptableRange
### Title: Acceptable Range Data
### Aliases: AcceptableRange
### Keywords: datasets

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("AcceptableRange", package = "DspikeIn")
  head(AcceptableRange)
  summary(AcceptableRange$Percentage)
}



cleanEx()
nameEx("MG_shapes")
### * MG_shapes

flush(stderr()); flush(stdout())

### Name: MG_shapes
### Title: Predefined Shape Vector for Plot Styling
### Aliases: MG_shapes
### Keywords: datasets

### ** Examples

# Example usage of MG_shapes in ggplot2
ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, shape = factor(cyl))) +
  ggplot2::geom_point(size = 4) +
  ggplot2::scale_shape_manual(values = MG_shapes[seq_len(3)]) +
  my_custom_theme()




cleanEx()
nameEx("Pre_processing_hashcodes")
### * Pre_processing_hashcodes

flush(stderr()); flush(stdout())

### Name: Pre_processing_hashcodes
### Title: Pre-process phyloseq or TSE object based on hashcodes
### Aliases: Pre_processing_hashcodes

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Subset phyloseq object to include only Tetragenococcus species
  Tetragenococcus <- phyloseq::subset_taxa(
    physeq_16SOTU,
    Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp."
  )

  # Extract hashcodes (OTU IDs)
  hashcodes <- row.names(phyloseq::otu_table(Tetragenococcus))

  # Temp file paths
  sum_prefix <- file.path(tempdir(), "merged_physeq_sum")
  max_prefix <- file.path(tempdir(), "merged_physeq_max")

  # Process data using the "sum" merging method
  processed_data_sum <- Pre_processing_hashcodes(
    physeq_16SOTU,
    hashcodes,
    merge_method = "sum",
    output_prefix = sum_prefix
  )

  # Convert to TreeSummarizedExperiment (TSE)
  Tetragenococcus_TSE <- convert_phyloseq_to_tse(Tetragenococcus)

  # Extract hashcodes for TSE
  hashcodes <- rownames(get_otu_table(Tetragenococcus_TSE))

  # Convert full phyloseq object to TSE format
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Load precomputed TSE object (optional)
  data("tse_16SOTU", package = "DspikeIn")

  # Process data using the "max" merging method
  processed_data_max <- Pre_processing_hashcodes(
    tse_16SOTU,
    hashcodes,
    merge_method = "max",
    output_prefix = max_prefix
  )

  # Clean up
  files_to_delete <- paste0(c(sum_prefix, max_prefix), "_processed.rds")
  invisible(file.remove(files_to_delete[file.exists(files_to_delete)]))
}



cleanEx()
nameEx("Pre_processing_species")
### * Pre_processing_species

flush(stderr()); flush(stdout())

### Name: Pre_processing_species
### Title: Pre-process taxa in a phyloseq or TSE object by merging
###   ASVs/OTUs
### Aliases: Pre_processing_species

### ** Examples

library(DspikeIn)
data("physeq_16SOTU", package = "DspikeIn")

species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")

# Merge species in phyloseq format
merged_sum <- Pre_processing_species(
  physeq_16SOTU,
  species_name,
  merge_method = "sum"
)

# Convert phyloseq to TSE format
tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

# Merge species in TSE format and write to tempdir
output_rds <- file.path(tempdir(), "merged_TSE_sum.rds")

merged_TSE_sum <- Pre_processing_species(
  tse_16SOTU,
  species_name,
  merge_method = "sum",
  output_file = output_rds
)



cleanEx()
nameEx("Pre_processing_species_list")
### * Pre_processing_species_list

flush(stderr()); flush(stdout())

### Name: Pre_processing_species_list
### Title: Preprocess and Merge Spike-in Species in a Phyloseq or TSE
###   Object
### Aliases: Pre_processing_species_list

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("tse", package = "DspikeIn")
  data("physeq", package = "DspikeIn")

  spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")

  merged_TSE <- Pre_processing_species_list(
    tse,
    spiked_species = spiked_species,
    merge_method = "sum"
  )

  merged_physeq <- Pre_processing_species_list(
    physeq,
    spiked_species = spiked_species,
    merge_method = "sum"
  )
}




cleanEx()
nameEx("RandomForest_selected")
### * RandomForest_selected

flush(stderr()); flush(stdout())

### Name: RandomForest_selected
### Title: Select Important ASVs/OTUs Using Random Forest
### Aliases: RandomForest_selected

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Perform Random Forest feature selection
  rf_physeq <- RandomForest_selected(
    physeq_16SOTU,
    response_var = "Host.genus",
    na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet")
  )

  # Load TreeSummarizedExperiment (TSE) object
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Perform Random Forest feature selection on TSE object
  rf_tse <- RandomForest_selected(
    tse_16SOTU,
    response_var = "Host.genus",
    na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet")
  )
}




cleanEx()
nameEx("adjust_abundance_one_third")
### * adjust_abundance_one_third

flush(stderr()); flush(stdout())

### Name: adjust_abundance_one_third
### Title: Adjust Abundance by a Factor
### Aliases: adjust_abundance_one_third

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")
  adjusted_physeq <- adjust_abundance_one_third(physeq_16SOTU, factor = 3)

  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  adjusted_tse <- adjust_abundance_one_third(tse_16SOTU, factor = 3)
}




cleanEx()
nameEx("adjusted_prevalence")
### * adjusted_prevalence

flush(stderr()); flush(stdout())

### Name: adjusted_prevalence
### Title: Adjust Prevalence in a Microbiome Object
### Aliases: adjusted_prevalence

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  ## Adjust prevalence in phyloseq
  adjusted_physeq <- adjusted_prevalence(physeq_16SOTU, method = "mean")

  ## Convert to TreeSummarizedExperiment
  tse_obj <- convert_phyloseq_to_tse(physeq_16SOTU)

  ## Adjust prevalence in TSE
  adjusted_tse <- adjusted_prevalence(tse_obj, method = "median")
}




cleanEx()
nameEx("alluvial_plot")
### * alluvial_plot

flush(stderr()); flush(stdout())

### Name: alluvial_plot
### Title: Generate an Alluvial Plot for Microbiome Data
### Aliases: alluvial_plot

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Convert phyloseq object to long format
  pps_Abs <- get_long_format_data(physeq_16SOTU)

  # Example of total reads calculation for relative abundance
  total_reads <- sum(pps_Abs$Abundance)
  print(paste("Total reads:", total_reads))

  # Generate an alluvial plot for relative abundance using the extended palette
  alluvial_plot_rel <- alluvial_plot(
    data = pps_Abs,
    axes = c("Env.broad.scale", "Host.genus", "Diet"),
    abundance_threshold = 0.01,
    fill_variable = "Phylum",
    silent = TRUE,
    abundance_type = "relative",
    top_taxa = 10,
    text_size = 4,
    legend_ncol = 1,
    custom_colors = DspikeIn::color_palette$cool_MG # Extended color palette
  )

  print(alluvial_plot_rel)

  # Convert phyloseq object to TreeSummarizedExperiment (TSE)
  tse_data <- convert_phyloseq_to_tse(physeq_16SOTU)
  tse_long <- get_long_format_data(tse_data)

  # Generate an alluvial plot for absolute abundance
  alluvial_plot_abs <- alluvial_plot(
    data = tse_long,
    axes = c("Env.broad.scale", "Host.genus", "Diet"),
    abundance_threshold = 8000,
    fill_variable = "Phylum",
    silent = TRUE,
    abundance_type = "absolute",
    top_taxa = 10,
    text_size = 4,
    legend_ncol = 1,
    custom_colors = DspikeIn::color_palette$cool_MG
  )

  print(alluvial_plot_abs)
}




cleanEx()
nameEx("calculate_list_average_scaling_factors")
### * calculate_list_average_scaling_factors

flush(stderr()); flush(stdout())

### Name: calculate_list_average_scaling_factors
### Title: Calculate Sample-specific Average Scaling Factors for Multiple
###   Spike-in Groups
### Aliases: calculate_list_average_scaling_factors

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq", package = "DspikeIn")

  # Step 1: Define the spike-in species groups and associated cell counts
  spiked_species_list <- list(
    c("Pseudomonas aeruginosa"),
    c("Escherichia coli"),
    c("Clostridium difficile")
  )

  spiked_cells_list <- c(10000, 20000, 15000)

  # Step 2: Compute the scaling factors
  scaling_factors <- calculate_list_average_scaling_factors(
    physeq,
    spiked_species_list,
    spiked_cells_list,
    merge_method = "sum"
  )

  # Step 3: Inspect scaling factors
  print(scaling_factors)
}



cleanEx()
nameEx("calculate_spikeIn_factors")
### * calculate_spikeIn_factors

flush(stderr()); flush(stdout())

### Name: calculate_spikeIn_factors
### Title: Calculate Scaling Factors for Spiked Species in Phyloseq or TSE
###   Object
### Aliases: calculate_spikeIn_factors

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE) &&
  requireNamespace("phyloseq", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  spiked_cells <- 1847
  species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
  merged_spiked_species <- "Tetragenococcus_halophilus"

  # --- Phyloseq Example ---
  spiked_16S_OTU <- phyloseq::subset_samples(physeq_16SOTU, spiked.volume %in% c("2", "1"))
  temp_output_file <- file.path(tempdir(), "merged_physeq_sum.rds")
  output_dir <- file.path(tempdir(), "spikeIn_factors_output")

  merged_physeq_sum <- Pre_processing_species(
    spiked_16S_OTU,
    species_name = merged_spiked_species,
    merge_method = "sum",
    output_file = temp_output_file
  )

  result_physeq <- calculate_spikeIn_factors(
    merged_physeq_sum,
    spiked_cells = spiked_cells,
    merged_spiked_species = merged_spiked_species,
    output_path = output_dir
  )

  print(result_physeq$scaling_factors)

  if (file.exists(temp_output_file)) unlink(temp_output_file, force = TRUE)
  if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)

  # --- TSE Example ---
  tse_data <- convert_phyloseq_to_tse(physeq_16SOTU)
  merged_tse_sum <- Pre_processing_species(
    tse_data,
    species_name = merged_spiked_species,
    merge_method = "sum"
  )

  result_tse <- calculate_spikeIn_factors(
    merged_tse_sum,
    spiked_cells = spiked_cells,
    merged_spiked_species = merged_spiked_species
  )

  print(result_tse$scaling_factors)
}



cleanEx()
nameEx("calculate_spike_percentage")
### * calculate_spike_percentage

flush(stderr()); flush(stdout())

### Name: calculate_spike_percentage
### Title: Calculate Spike Percentage for Specified Taxa in a Phyloseq or
###   TSE Object
### Aliases: calculate_spike_percentage

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  # Load example phyloseq object
  data("physeq_16SOTU", package = "DspikeIn")

  # ----------- Phyloseq Example -----------
  species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
  merged_spiked_species <- "Tetragenococcus_halophilus"

  # Pre-process the phyloseq object to merge spike-in taxa
  merged_physeq <- Pre_processing_species(
    physeq_16SOTU,
    species_name = species_name,
    merge_method = "sum"
  )

  # Perform spike-in percentage calculation
  output_docx <- file.path(tempdir(), "spike_summary_physeq.docx")
  result_physeq <- calculate_spike_percentage(
    obj = merged_physeq,
    merged_spiked_species = merged_spiked_species,
    output_file = output_docx,
    passed_range = c(0.1, 20)
  )
  print(result_physeq)

  # ----------- TreeSummarizedExperiment Example -----------
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  merged_tse <- Pre_processing_species(
    tse_16SOTU,
    species_name = species_name,
    merge_method = "sum"
  )

  output_docx_tse <- file.path(tempdir(), "spike_summary_tse.docx")
  result_tse <- calculate_spike_percentage(
    obj = merged_tse,
    merged_spiked_species = merged_spiked_species,
    output_file = output_docx_tse,
    passed_range = c(0.1, 20)
  )
  print(result_tse)

  # Clean up temporary files
  if (file.exists(output_docx)) unlink(output_docx, force = TRUE)
  if (file.exists(output_docx_tse)) unlink(output_docx_tse, force = TRUE)
}



cleanEx()
nameEx("calculate_spike_percentage_list")
### * calculate_spike_percentage_list

flush(stderr()); flush(stdout())

### Name: calculate_spike_percentage_list
### Title: Calculate Spike-in Percentage for Specified Taxa
### Aliases: calculate_spike_percentage_list

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq", package = "DspikeIn")
  spiked_species_list <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")

  temp_docx <- file.path(tempdir(), "merged_result.docx")
  temp_csv <- sub(".docx", ".csv", temp_docx)

  result <- calculate_spike_percentage_list(
    obj = physeq,
    merged_spiked_species = spiked_species_list,
    output_path = temp_docx,
    passed_range = c(0.1, 10)
  )

  print(result)

  # Clean up
  if (file.exists(temp_docx)) unlink(temp_docx, force = TRUE)
  if (file.exists(temp_csv)) unlink(temp_csv, force = TRUE)
}



cleanEx()
nameEx("calculate_summary_stats_table")
### * calculate_summary_stats_table

flush(stderr()); flush(stdout())

### Name: calculate_summary_stats_table
### Title: Calculate Summary Statistics Table
### Aliases: calculate_summary_stats_table

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  # --- Phyloseq example ---
  data("physeq_16SOTU", package = "DspikeIn")
  absolute_count <- phyloseq::otu_table(physeq_16SOTU)

  tmp_docx <- file.path(tempdir(), "physeq_summary.docx")
  summary_table_physeq <- calculate_summary_stats_table(
    data = as.data.frame(absolute_count),
    output_path = tmp_docx
  )
  print(summary_table_physeq)
  if (file.exists(tmp_docx)) file.remove(tmp_docx)

  # --- TSE example ---
  data("tse", package = "DspikeIn")
  tse_counts <- SummarizedExperiment::assay(tse)

  tmp_docx2 <- file.path(tempdir(), "tse_summary.docx")
  summary_table_tse <- calculate_summary_stats_table(
    data = as.data.frame(tse_counts),
    output_path = tmp_docx2
  )
  print(summary_table_tse)
  if (file.exists(tmp_docx2)) file.remove(tmp_docx2)
}




cleanEx()
nameEx("color_palette")
### * color_palette

flush(stderr()); flush(stdout())

### Name: color_palette
### Title: Original, Extended, and Nature-Inspired Color Palette Sequence
### Aliases: color_palette
### Keywords: datasets

### ** Examples

# Example using the Mar_palette
ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, color = factor(cyl))) +
  ggplot2::geom_point(size = 4) +
  ggplot2::scale_color_manual(values = color_palette$Mar_palette) +
  ggplot2::theme_minimal()




cleanEx()
nameEx("conclusion")
### * conclusion

flush(stderr()); flush(stdout())

### Name: conclusion
### Title: Compute Summary Statistics for Spiked Species
### Aliases: conclusion

### ** Examples

## -------------------------------
## Example 1: Using phyloseq object
## -------------------------------
library(DspikeIn)
data("physeq_16SOTU", package = "DspikeIn")

# Merge spike-in species
species_name <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
merged_sum <- Pre_processing_species(
  obj = physeq_16SOTU,
  species_name = species_name,
  merge_method = "sum"
)

# Compute summary statistics
output_doc <- file.path(tempdir(), "summary_phyloseq.docx")

results_physeq <- conclusion(
  obj = merged_sum,
  merged_spiked_species = "Tetragenococcus_halophilus",
  max_passed_range = 20,
  output_path = output_doc
)
print(results_physeq$summary_stats)

## -----------------------------------------------
## Example 2: Using TreeSummarizedExperiment object
## -----------------------------------------------
tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

output_doc_tse <- file.path(tempdir(), "summary_tse.docx")
results_tse <- conclusion(
  obj = tse_16SOTU,
  merged_spiked_species = "Tetragenococcus_halophilus",
  max_passed_range = 20,
  output_path = output_doc_tse
)
print(results_tse$summary_stats)




cleanEx()
nameEx("convert_categorical_to_factors")
### * convert_categorical_to_factors

flush(stderr()); flush(stdout())

### Name: convert_categorical_to_factors
### Title: Convert Categorical Columns to Factors in Sample Data
### Aliases: convert_categorical_to_factors

### ** Examples

data("physeq_16SOTU", package = "DspikeIn")
ps_factor <- convert_categorical_to_factors(physeq_16SOTU)



cleanEx()
nameEx("convert_phyloseq_to_tse")
### * convert_phyloseq_to_tse

flush(stderr()); flush(stdout())

### Name: convert_phyloseq_to_tse
### Title: Convert a 'phyloseq' Object to a 'TreeSummarizedExperiment'
### Aliases: convert_phyloseq_to_tse

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Convert phyloseq object to TreeSummarizedExperiment
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
}



cleanEx()
nameEx("convert_to_absolute_counts")
### * convert_to_absolute_counts

flush(stderr()); flush(stdout())

### Name: convert_to_absolute_counts
### Title: Convert Relative ASV/OTU Counts to Absolute Counts
### Aliases: convert_to_absolute_counts

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  spiked_cells <- 1847
  species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp.")
  merged_spiked_species <- "Tetragenococcus_halophilus"

  spiked_16S_OTU <- phyloseq::subset_samples(physeq_16SOTU, spiked.volume %in% c("2", "1"))
  Spiked_16S_sum_scaled <- Pre_processing_species(
    spiked_16S_OTU,
    species_name,
    merge_method = "sum",
    output_file = file.path(tempdir(), "merged_physeq_sum.rds")
  )

  result <- calculate_spikeIn_factors(
    Spiked_16S_sum_scaled,
    spiked_cells,
    merged_spiked_species
  )
  scaling_factors <- result$scaling_factors

  result_physeq <- convert_to_absolute_counts(
    Spiked_16S_sum_scaled,
    scaling_factors,
    output_dir = tempdir()
  )
  abs_counts_physeq <- result_physeq$absolute_counts
  physeq_adj <- result_physeq$obj_adj

  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  spiked_16S_OTU_TSE <- tse_16SOTU[, tse_16SOTU$spiked.volume %in% c("2", "1")]

  Spiked_16S_sum_tse_scaled <- Pre_processing_species(
    spiked_16S_OTU_TSE,
    species_name,
    merge_method = "sum",
    output_file = file.path(tempdir(), "merged_tse_sum.rds")
  )

  result <- calculate_spikeIn_factors(
    Spiked_16S_sum_tse_scaled,
    spiked_cells,
    merged_spiked_species
  )
  scaling_factors <- result$scaling_factors

  result_tse <- convert_to_absolute_counts(
    spiked_16S_OTU_TSE,
    scaling_factors,
    output_dir = tempdir()
  )
  abs_counts_tse <- result_tse$absolute_counts
  tse_adj <- result_tse$obj_adj
}



cleanEx()
nameEx("convert_tse_to_phyloseq")
### * convert_tse_to_phyloseq

flush(stderr()); flush(stdout())

### Name: convert_tse_to_phyloseq
### Title: Convert a 'TreeSummarizedExperiment' to a 'phyloseq' Object
### Aliases: convert_tse_to_phyloseq

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Convert phyloseq to TSE
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Convert TSE back to phyloseq
  phy_M <- convert_tse_to_phyloseq(tse_16SOTU)
  print(phy_M)
}



cleanEx()
nameEx("create_directory")
### * create_directory

flush(stderr()); flush(stdout())

### Name: create_directory
### Title: Create a Directory and Optionally Set as Working Directory
### Aliases: create_directory

### ** Examples

if (interactive()) {
  # Save the current working directory
  old_wd <- getwd()

  # Use a temporary directory for safe example use
  tmp_dir <- file.path(tempdir(), "example_new_dir")

  # Create the directory and set it as the working directory
  create_directory(tmp_dir, set_working_dir = TRUE)

  # Do something inside tmp_dir...

  # Restore the original working directory
  setwd(old_wd)

  # Remove the created directory
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
}



cleanEx()
nameEx("degree_network")
### * degree_network

flush(stderr()); flush(stdout())

### Name: degree_network
### Title: Analyze and Visualize a Microbial Network
### Aliases: degree_network

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  Complete <- load_graphml(system.file("extdata", "Complete.graphml", package = "DspikeIn"))

  # Save metrics to a temporary file
  temp_metrics <- file.path(tempdir(), "Global_Network_Metrics.csv")

  result <- degree_network(
    graph_path = Complete,
    save_metrics = TRUE,
    metrics_path = temp_metrics
  )

  print(result$metrics)
  print(result$plot)

  # Clean up temporary file
  unlink(temp_metrics)
}




cleanEx()
nameEx("detect_common_asvs_taxa")
### * detect_common_asvs_taxa

flush(stderr()); flush(stdout())

### Name: detect_common_asvs_taxa
### Title: Detect Common ASVs and Taxa from Multiple Phyloseq or TSE
###   Objects
### Aliases: detect_common_asvs_taxa

### ** Examples

## Not run: 
##D if (requireNamespace("DspikeIn", quietly = TRUE)) {
##D   # Example with phyloseq objects
##D   common_physeq <- detect_common_asvs_taxa(list(physeq1, physeq2, physeq3))
##D 
##D   # Example with TreeSummarizedExperiment objects
##D   common_tse <- detect_common_asvs_taxa(list(tse1, tse2, tse3))
##D }
## End(Not run)



cleanEx()
nameEx("extract_neighbors")
### * extract_neighbors

flush(stderr()); flush(stdout())

### Name: extract_neighbors
### Title: Extract First and Second Neighbors of a Target Node
### Aliases: extract_neighbors

### ** Examples

# Load the built-in Complete graph
complete_graph <- load_graphml("Complete.graphml")
result1 <- extract_neighbors(
  graph = complete_graph,
  target_node = "OTU69:Basidiobolus_sp"
)
print(result1$summary)


# Load from an external GraphML file (ensure the file path is correct)
# external_graph <- load_graphml("~/custom_network.graphml")
# extract_neighbors(external_graph, target_node = "SomeNode")




cleanEx()
nameEx("filter_and_split_abundance")
### * filter_and_split_abundance

flush(stderr()); flush(stdout())

### Name: filter_and_split_abundance
### Title: Filter and Split Abundance Data by Threshold
### Aliases: filter_and_split_abundance

### ** Examples

data("physeq_ITSOTU", package = "DspikeIn")

# Return results without saving
output <- filter_and_split_abundance(physeq_ITSOTU, threshold = 0.05)

# With saving
tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)

output <- filter_and_split_abundance(tse_ITSOTU,
  threshold = 0.05,
  output_prefix = file.path(tempdir(), "abund")
)




cleanEx()
nameEx("get_long_format_data")
### * get_long_format_data

flush(stderr()); flush(stdout())

### Name: get_long_format_data
### Title: Convert a Phyloseq or TSE Object into a Long-Format Data Frame
### Aliases: get_long_format_data

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Turn it to the long format
  melted_ph <- get_long_format_data(physeq_16SOTU)
  melted <- get_long_format_data(tse_16SOTU)
}



cleanEx()
nameEx("gm_mean")
### * gm_mean

flush(stderr()); flush(stdout())

### Name: gm_mean
### Title: Calculate Geometric Mean
### Aliases: gm_mean
### Keywords: -----------------------------------------------------------
###   NULL microbiome normalization

### ** Examples

vec <- c(1, 10, 100, 1000)
gm_mean(vec)



cleanEx()
nameEx("metadata_full")
### * metadata_full

flush(stderr()); flush(stdout())

### Name: metadata_full
### Title: Metadata for Microbiome Samples
### Aliases: metadata_full
### Keywords: datasets

### ** Examples

data(metadata_full)
head(metadata_full)
summary(metadata_full)



cleanEx()
nameEx("my_custom_theme")
### * my_custom_theme

flush(stderr()); flush(stdout())

### Name: my_custom_theme
### Title: Custom ggplot2 Theme with Consistent Aesthetics
### Aliases: my_custom_theme

### ** Examples

library(ggplot2)
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(size = 3) +
  my_custom_theme()
print(p)




cleanEx()
nameEx("node_level_metrics")
### * node_level_metrics

flush(stderr()); flush(stdout())

### Name: node_level_metrics
### Title: Compute and Visualize Node-Level Network Metrics
### Aliases: node_level_metrics

### ** Examples

library(igraph)
set.seed(42)
# For external graphml please use full address
# Load internal graphml
Complete <- load_graphml("Complete.graphml")

# Compute node-level metrics
result <- node_level_metrics(Complete)

# View computed metrics
print(result$metrics)

# Show the first 4x4 plot
print(result$plots$plot1)

# Show the second 4x4 plot
print(result$plots$plot2)

# Show facet plot
print(result$facet_plot)

# Print metrics and flextable
print(result$metrics)
print(result$flextable)



cleanEx()
nameEx("norm.DESeq")
### * norm.DESeq

flush(stderr()); flush(stdout())

### Name: norm.DESeq
### Title: DESeq Normalization with Pseudocount and Integer Conversion
### Aliases: norm.DESeq

### ** Examples




cleanEx()
nameEx("normalization_set")
### * normalization_set

flush(stderr()); flush(stdout())

### Name: normalization_set
### Title: Apply the Selected Normalization Method to the Phyloseq and TSE
###   Objects
### Aliases: normalization_set

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU

  # TC Normalization
  result_TC <- normalization_set(ps, method = "TC", groups = "Host.species")
  normalized_ps_TC <- result_TC$dat.normed
  scaling_factors_TC <- result_TC$scaling.factor

  # UQ Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_UQ <- normalization_set(ps, method = "UQ", groups = "Host.species")
  normalized_ps_UQ <- result_UQ$dat.normed
  scaling_factors_UQ <- result_UQ$scaling.factor

  # Median Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_med <- normalization_set(ps, method = "med", groups = "Host.species")
  normalized_ps_med <- result_med$dat.normed
  scaling_factors_med <- result_med$scaling.factor

  # DESeq Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  ps_n <- remove_zero_negative_count_samples(ps)
  result_DESeq <- normalization_set(ps_n, method = "DESeq", groups = "Animal.type")
  normalized_ps_DESeq <- result_DESeq$dat.normed
  scaling_factors_DESeq <- result_DESeq$scaling.factor

  # Poisson Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_Poisson <- normalization_set(ps, method = "Poisson", groups = "Host.genus")
  normalized_ps_Poisson <- result_Poisson$dat.normed
  scaling_factors_Poisson <- result_Poisson$scaling.factor

  # Quantile Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_QN <- normalization_set(ps, method = "QN")
  normalized_ps_QN <- result_QN$dat.normed
  scaling_factors_QN <- result_QN$scaling.factor

  # TMM Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_TMM <- normalization_set(ps, method = "TMM", groups = "Animal.type")
  normalized_ps_TMM <- result_TMM$dat.normed
  scaling_factors_TMM <- result_TMM$scaling.factor

  # CLR Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_clr <- normalization_set(ps, method = "clr")
  normalized_ps_clr <- result_clr$dat.normed
  scaling_factors_clr <- result_clr$scaling.factor

  # Rarefying
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_rar <- normalization_set(ps, method = "rar")
  normalized_ps_rar <- result_rar$dat.normed
  scaling_factors_rar <- result_rar$scaling.factor

  # CSS Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_css <- normalization_set(ps, method = "css")
  normalized_ps_css <- result_css$dat.normed
  scaling_factors_css <- result_css$scaling.factor

  # TSS Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_tss <- normalization_set(ps, method = "tss")
  normalized_ps_tss <- result_tss$dat.normed
  scaling_factors_tss <- result_tss$scaling.factor

  # RLE Normalization
  data("physeq_16SOTU", package = "DspikeIn")
  ps <- physeq_16SOTU
  result_rle <- normalization_set(ps, method = "rle")
  normalized_ps_rle <- result_rle$dat.normed
  scaling_factors_rle <- result_rle$scaling.factor
}




cleanEx()
nameEx("perform_and_visualize_DA")
### * perform_and_visualize_DA

flush(stderr()); flush(stdout())

### Name: perform_and_visualize_DA
### Title: Perform and Visualize Differential Abundance Analysis with edgeR
###   or DESeq2
### Aliases: perform_and_visualize_DA

### ** Examples

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Single contrast (DESeq2)
  res_single <- perform_and_visualize_DA(
    obj = physeq_16SOTU,
    method = "DESeq2",
    group_var = "Diet",
    contrast = c("Insectivore", "Carnivore"),
    target_glom = "Genus"
  )

  # Multiple contrasts for one group_var
  contrast_list <- list(
    c("Insectivore", "Carnivore"),
    c("Omnivore", "Herbivore")
  )

  res_multi <- perform_and_visualize_DA(
    obj = physeq_16SOTU,
    method = "DESeq2",
    group_var = "Diet",
    contrast = contrast_list,
    global_fdr = TRUE
  )

  # Named contrast list for multiple group_vars
  contrast_named <- list(
    Diet = list(
      c("Insectivore", "Carnivore"),
      c("Omnivore", "Carnivore")
    ),
    Animal.type = list(
      c("Frog", "Salamander")
    )
  )

  res_multi_factor <- perform_and_visualize_DA(
    obj = physeq_16SOTU,
    method = "DESeq2",
    significance_level = 0.01,
    contrast = contrast_named,
    target_glom = "Genus",
    global_fdr = TRUE
  )

  # Create a combined factor of Animal.type and Diet
  phyloseq::sample_data(physeq_16SOTU)$ComboGroup <- factor(interaction(
    phyloseq::sample_data(physeq_16SOTU)$Animal.type,
    phyloseq::sample_data(physeq_16SOTU)$Diet,
    drop = TRUE
  ))

  # Define valid contrasts with sufficient sample sizes
  contrast_list <- list(
    c("Salamander.Insectivore", "Lizard.Insectivore"),
    c("Salamander.Carnivore", "Snake.Carnivore"),
    c("Salamander.Carnivore", "Frog.Carnivore")
  )

  # Perform multi-contrast DA analysis on the ComboGroup
  res_combo <- perform_and_visualize_DA(
    obj = physeq_16SOTU,
    method = "DESeq2",
    group_var = "ComboGroup",
    contrast = contrast_list,
    target_glom = "Genus",
    global_fdr = TRUE
  )

  # Example: view one of the plots
  res_combo$Frog.Insectivore_vs_Lizard.Insectivore$bar_plot
}




cleanEx()
nameEx("physeq")
### * physeq

flush(stderr()); flush(stdout())

### Name: physeq
### Title: General Example Phyloseq Object
### Aliases: physeq
### Keywords: datasets

### ** Examples

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(physeq)
  physeq
  phyloseq::sample_names(physeq)
  phyloseq::taxa_names(physeq)
  summary(physeq)
}



cleanEx()
nameEx("physeq_16SOTU")
### * physeq_16SOTU

flush(stderr()); flush(stdout())

### Name: physeq_16SOTU
### Title: Example Phyloseq Object for 16S OTUs
### Aliases: physeq_16SOTU
### Keywords: datasets

### ** Examples

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(physeq_16SOTU)
  physeq_16SOTU
  summary(physeq_16SOTU)
  phyloseq::sample_names(physeq_16SOTU)
  phyloseq::taxa_names(physeq_16SOTU)
}



cleanEx()
nameEx("physeq_ITSOTU")
### * physeq_ITSOTU

flush(stderr()); flush(stdout())

### Name: physeq_ITSOTU
### Title: Example Phyloseq Object for ITS OTUs
### Aliases: physeq_ITSOTU
### Keywords: datasets

### ** Examples

if (requireNamespace("phyloseq", quietly = TRUE)) {
  data(physeq_ITSOTU)
  physeq_ITSOTU
  summary(physeq_ITSOTU)
  phyloseq::sample_data(physeq_ITSOTU)
  phyloseq::taxa_names(physeq_ITSOTU)
}



cleanEx()
nameEx("plot_core_microbiome_custom")
### * plot_core_microbiome_custom

flush(stderr()); flush(stdout())

### Name: plot_core_microbiome_custom
### Title: Plot Core Microbiome Prevalence Heatmap (Phyloseq & TSE
###   Compatible)
### Aliases: plot_core_microbiome_custom

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Define custom detection parameters
  custom_detections <- list(
    prevalences = seq(0.03, 1, 0.01),
    thresholds = 10^seq(log10(0.03), log10(1), length = 10),
    min_prevalence = 0.3,
    taxa_order = "ascending"
  )

  # Define output paths using tempdir()
  core_csv <- file.path(tempdir(), "core_microbiome.csv")
  core_rds <- file.path(tempdir(), "core_microbiome.rds")

  # Generate a core microbiome plot with custom settings
  plot_result <- plot_core_microbiome_custom(
    obj = physeq_16SOTU,
    detections = custom_detections,
    taxrank = "Genus",
    output_core_rds = core_rds,
    output_core_csv = core_csv
  )

  # Print the resulting plot
  print(plot_result)
}



cleanEx()
nameEx("plot_spikein_tree_diagnostic")
### * plot_spikein_tree_diagnostic

flush(stderr()); flush(stdout())

### Name: plot_spikein_tree_diagnostic
### Title: Spike-in Tree Diagnostic Plot
### Aliases: plot_spikein_tree_diagnostic

### ** Examples

## Not run: 
##D if (
##D   requireNamespace("DspikeIn", quietly = TRUE) &&
##D     requireNamespace("phyloseq", quietly = TRUE) &&
##D     requireNamespace("TreeSummarizedExperiment", quietly = TRUE) &&
##D     requireNamespace("ggplot2", quietly = TRUE) &&
##D     requireNamespace("ggtree", quietly = TRUE) &&
##D     requireNamespace("ggtreeExtra", quietly = TRUE) &&
##D     requireNamespace("ggnewscale", quietly = TRUE)
##D ) {
##D   # Load synthetic test dataset from DspikeIn
##D   data("physeq_16SOTU", package = "DspikeIn")
##D 
##D   # Filter to a known spike-in genus
##D   spikein <- phyloseq::subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
##D 
##D   # Plot diagnostic using phyloseq object
##D   plot_spikein_tree_diagnostic(
##D     obj = spikein,
##D     metadata_var = "Animal.type",
##D     save_plot = FALSE
##D   )
##D 
##D   # Convert to TreeSummarizedExperiment object
##D   tse_spikein <- convert_phyloseq_to_tse(spikein)
##D 
##D   # Plot diagnostic using TSE object
##D   plot_spikein_tree_diagnostic(
##D     obj = tse_spikein,
##D     metadata_var = "Animal.type",
##D     save_plot = FALSE
##D   )
##D }
## End(Not run)



cleanEx()
nameEx("plotbar_abundance")
### * plotbar_abundance

flush(stderr()); flush(stdout())

### Name: plotbar_abundance
### Title: Taxa Bar Plot Without Aggregation (Relative or Absolute
###   Abundance)
### Aliases: plotbar_abundance

### ** Examples

# Load required package
if (requireNamespace("phyloseq", quietly = TRUE)) {

  # Load example data
  data("physeq_ITSOTU", package = "DspikeIn")

  # Subset: Eurycea salamanders from Blue Ridge, exclude unwanted genera
  Des <- physeq_ITSOTU |>
    phyloseq::subset_taxa(Genus != "Dekkera") |>
    phyloseq::subset_samples(Clade.Order == "Caudate") |>
    phyloseq::subset_samples(Host.genus == "Eurycea") |>
    phyloseq::subset_samples(Ecoregion.III == "Blue Ridge")

  # Clean taxa: remove NAs or blanks in Phylum, filter low-abundance
  Des_filtered <- phyloseq::subset_taxa(Des, !is.na(Phylum) & Phylum != "")
  Des_ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(Des_filtered) > 99, Des_filtered)

  # Plot taxa abundance with full control
  plotbar_abundance(
    physeq = Des_ps,
    normalize = TRUE,
    treatment_variable = "Diet",
    abundance_type = "absolute",
    x_angle = 0,
    fill_variable = "Phylum",
    palette = DspikeIn::color_palette$mix_MG,
    legend_size = 10,
    legend_columns = 1,
    x_scale = "free",
    xlab = NULL
  )
}




cleanEx()
nameEx("proportion_adj")
### * proportion_adj

flush(stderr()); flush(stdout())

### Name: proportion_adj
### Title: Proportionally Adjust Abundance
### Aliases: proportion_adj

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  # Load phyloseq object
  data("physeq_16SOTU", package = "DspikeIn")

  normalized_physeq <- proportion_adj(
    physeq_16SOTU,
    output_file = file.path(tempdir(), "proportion_adjusted_physeq.rds")
  )
  print(normalized_physeq)

  # Convert to TSE and apply
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  normalized_tse <- proportion_adj(
    tse_16SOTU,
    output_file = file.path(tempdir(), "proportion_adjusted_tse.rds")
  )
  print(normalized_tse)
}



cleanEx()
nameEx("quadrant_plot")
### * quadrant_plot

flush(stderr()); flush(stdout())

### Name: quadrant_plot
### Title: Generate Custom Quadrant Plots for Node Metrics
### Aliases: quadrant_plot

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  g <- load_graphml("Complete.graphml")

  # Compute node-level metrics
  result <- node_level_metrics(g)
  metrics <- result$metrics

  # Generate a quadrant plot using Degree and Efficiency
  plot <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Efficiency")
  print(plot)
}




cleanEx()
nameEx("random_subsample_WithReductionFactor")
### * random_subsample_WithReductionFactor

flush(stderr()); flush(stdout())

### Name: random_subsample_WithReductionFactor
### Title: Random Subsampling with Reduction Factor
### Aliases: random_subsample_WithReductionFactor

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE) &&
  requireNamespace("phyloseq", quietly = TRUE) &&
  requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")
  red <- random_subsample_WithReductionFactor(physeq_16SOTU, reduction_factor = 10)
  summary_stats <- summ_phyloseq_sampleID(red)
  print(summary_stats)

  tse <- convert_phyloseq_to_tse(physeq_16SOTU)
  red_tse <- random_subsample_WithReductionFactor(tse, reduction_factor = 10)
  summary_stats <- summ_phyloseq_sampleID(red_tse)
  print(summary_stats)
}




cleanEx()
nameEx("randomsubsample_Trimmed_evenDepth")
### * randomsubsample_Trimmed_evenDepth

flush(stderr()); flush(stdout())

### Name: randomsubsample_Trimmed_evenDepth
### Title: Subsampling to an Equal Sequencing Depth
### Aliases: randomsubsample_Trimmed_evenDepth

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_ITSOTU", package = "DspikeIn")
  tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
  rarefied <- randomsubsample_Trimmed_evenDepth(tse_ITSOTU, smalltrim = 0.001)
  print(rarefied)
}




cleanEx()
nameEx("regression_plot")
### * regression_plot

flush(stderr()); flush(stdout())

### Name: regression_plot
### Title: Create a Regression Plot with Faceting by Range
### Aliases: regression_plot

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("metadata_full", package = "DspikeIn")

  plot_object <- regression_plot(
    data = metadata_full,
    x_var = "Observed",
    y_var = "Total_Reads_spiked",
    custom_range = c(0.1, 15, 30, 50, 75, 100)
  )

  # Print the plot output
  print(plot_object)
}



cleanEx()
nameEx("relativized_filtered_taxa")
### * relativized_filtered_taxa

flush(stderr()); flush(stdout())

### Name: relativized_filtered_taxa
### Title: Filter Taxa from a Phyloseq or TSE Object Based on Custom
###   Thresholds
### Aliases: relativized_filtered_taxa

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Apply relative filtering on taxa
  FT <- relativized_filtered_taxa(
    physeq_16SOTU,
    threshold_percentage = 0.001,
    threshold_mean_abundance = 1,
    threshold_count = 5,
    threshold_relative_abundance = 0.001
  )
}



cleanEx()
nameEx("remove_zero_negative_count_samples")
### * remove_zero_negative_count_samples

flush(stderr()); flush(stdout())

### Name: remove_zero_negative_count_samples
### Title: Remove Samples with Zero, Negative Counts, or NA Values and Add
###   Pseudocount
### Aliases: remove_zero_negative_count_samples

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  library(DspikeIn)
  data("physeq_16SOTU", package = "DspikeIn")

  # Remove samples with zero/negative/NA counts and add pseudocount
  cleaned_ps <- remove_zero_negative_count_samples(
    physeq_16SOTU,
    pseudocount = 1e-6
  )
}



cleanEx()
nameEx("ridge_plot_it")
### * ridge_plot_it

flush(stderr()); flush(stdout())

### Name: ridge_plot_it
### Title: Generate Ridge Plots for Taxonomic Abundance Distribution
### Aliases: ridge_plot_it

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  # Load phyloseq object
  data("physeq_16SOTU", package = "DspikeIn")
  ridge_physeq <- ridge_plot_it(physeq_16SOTU, taxrank = "Family", top_n = 10)


  # convert phyloseq object to TSE
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  ridge_tse <- ridge_plot_it(tse_16SOTU, taxrank = "Family", top_n = 10)
}




cleanEx()
nameEx("simulate_network_robustness")
### * simulate_network_robustness

flush(stderr()); flush(stdout())

### Name: simulate_network_robustness
### Title: Simulate Network Robustness under Node Removal
### Aliases: simulate_network_robustness

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  Complete <- load_graphml("Complete.graphml")

  # Simulate robustness by removing 200 highest-degree nodes
  robustness_degree <- simulate_network_robustness(
    graph = Complete,
    steps = 200,
    removal_strategy = "degree"
  )

  # Simulate robustness with random node removal
  robustness_random <- simulate_network_robustness(
    graph = Complete,
    steps = 200,
    removal_strategy = "random"
  )

  # Simulate robustness with betweenness-based node removal
  robustness_betweenness <- simulate_network_robustness(
    graph = Complete,
    steps = 200,
    removal_strategy = "betweenness"
  )

  # Print robustness plots
  print(robustness_degree$plot)
  print(robustness_random$plot)
  print(robustness_betweenness$plot)
}



cleanEx()
nameEx("summ_ASV_OTUID")
### * summ_ASV_OTUID

flush(stderr()); flush(stdout())

### Name: summ_ASV_OTUID
### Title: Summarize ASV Data Based on ASV_ID
### Aliases: summ_ASV_OTUID

### ** Examples

# Example with a phyloseq object
if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_ITSOTU", package = "DspikeIn")
  summary_physeq <- summ_ASV_OTUID(physeq_ITSOTU)

  # Example with a TreeSummarizedExperiment object
  tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
  summary_tse <- summ_ASV_OTUID(tse_ITSOTU)
}




cleanEx()
nameEx("summ_count_phyloseq")
### * summ_count_phyloseq

flush(stderr()); flush(stdout())

### Name: summ_count_phyloseq
### Title: Summary Statistics of a Phyloseq or TSE Object
### Aliases: summ_count_phyloseq

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Summarize counts for the phyloseq object
  summary_stats_physeq <- summ_count_phyloseq(physeq_16SOTU)

  # Convert phyloseq object to a TreeSummarizedExperiment (TSE)
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)

  # Summarize counts for the TSE object
  summary_stats_tse <- summ_count_phyloseq(tse_16SOTU)
}




cleanEx()
nameEx("summ_phyloseq_sampleID")
### * summ_phyloseq_sampleID

flush(stderr()); flush(stdout())

### Name: summ_phyloseq_sampleID
### Title: Generate Summary Statistics for Each Sample
### Aliases: summ_phyloseq_sampleID

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  # Summarize the phyloseq object
  summary_stats_physeq <- summ_phyloseq_sampleID(physeq_16SOTU)
  print(summary_stats_physeq)

  # Convert to TreeSummarizedExperiment
  tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
  summary_stats_tse <- summ_phyloseq_sampleID(tse_16SOTU)
  print(summary_stats_tse)
}




cleanEx()
nameEx("taxa_barplot")
### * taxa_barplot

flush(stderr()); flush(stdout())

### Name: taxa_barplot
### Title: Generate a Taxa Barplot with Relative or Absolute Abundance
### Aliases: taxa_barplot

### ** Examples

# Example 1: Relative abundance barplot for Genus
data("physeq_16SOTU", package = "DspikeIn")
bp_rel <- taxa_barplot(
  physeq = physeq_16SOTU,
  target_glom = "Genus",
  fill_variable = "Family",
  treatment_variable = "Diet",
  abundance_type = "relative",
  top_n_taxa = 20,
  legend_size = 10,
  x_scale = "free",
  legend_columns = 1,
  palette = DspikeIn::color_palette$MG
)
print(bp_rel$barplot)

# Example 2: Absolute abundance barplot for Family with faceting
data("physeq_ITSOTU", package = "DspikeIn")
tse_ITSOTU <- convert_phyloseq_to_tse(physeq_ITSOTU)
bp_abs <- taxa_barplot(
  physeq = tse_ITSOTU,
  target_glom = "Genus",
  treatment_variable = "Animal.type",
  abundance_type = "absolute",
  fill_variable = "Family",
  facet_variable = "Diet",
  top_n_taxa = 10,
  x_scale = "fixed",
  xlab = NULL,
  legend_columns = 1,
  x_angle = 25,
  palette = DspikeIn::color_palette$cool_MG,
  legend_size = 12
)
print(bp_abs$barplot)




cleanEx()
nameEx("tidy_phyloseq_tse")
### * tidy_phyloseq_tse

flush(stderr()); flush(stdout())

### Name: tidy_phyloseq_tse
### Title: Tidy a Phyloseq or TreeSummarizedExperiment Object
### Aliases: tidy_phyloseq_tse

### ** Examples

if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")
  tidy_physeq <- tidy_phyloseq_tse(physeq_16SOTU)
}


if (requireNamespace("DspikeIn", quietly = TRUE)) {
  data("physeq_16SOTU", package = "DspikeIn")

  tidy_physeq <- tidy_phyloseq_tse(physeq_16SOTU)
}




cleanEx()
nameEx("tse")
### * tse

flush(stderr()); flush(stdout())

### Name: tse
### Title: Example TreeSummarizedExperiment (TSE) Object
### Aliases: tse
### Keywords: datasets

### ** Examples

if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
  data(tse)
  tse
  SummarizedExperiment::assay(tse)
  SummarizedExperiment::colData(tse)
  SummarizedExperiment::rowData(tse)
}



cleanEx()
nameEx("validate_spikein_clade")
### * validate_spikein_clade

flush(stderr()); flush(stdout())

### Name: validate_spikein_clade
### Title: Validate Spike-In Clade Consistency with NJ Tree and Bootstrap
### Aliases: validate_spikein_clade

### ** Examples

ref_fasta <- system.file("extdata", "Ref.fasta", package = "DspikeIn")
sample_fasta <- system.file("extdata", "Sample.fasta", package = "DspikeIn")
result <- validate_spikein_clade(ref_fasta, sample_fasta)




cleanEx()
nameEx("weight_Network")
### * weight_Network

flush(stderr()); flush(stdout())

### Name: weight_Network
### Title: Analyze and Visualize a Microbial Network
### Aliases: weight_Network

### ** Examples


if (requireNamespace("DspikeIn", quietly = TRUE)) {
  Complete <- load_graphml("Complete.graphml")

  # Run network weighting function on the loaded dataset
  # Load a specific GraphML dataset from DspikeIn
  result <- weight_Network(graph_path = "Complete.graphml")

  # Load a custom GraphML file from user directory,
  # for external graphml please use **full address**
  print(result$plot)

  # View network metrics
  result$metrics
  # Optional: Clean up generated metrics file if saved
  unlink("Global_Network_Metrics.csv")
}




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
