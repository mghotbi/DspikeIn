# DspikeIn

---

### DspikeIn Package
The **DspikeIn** package was developed to facilitate:
- Verifying the phylogenetic distances of ASVs/OTUs resulting from spiked species.
- Preprocessing data.
- Calculating the spike-in scaling factor.
- Converting relative abundance to absolute abundance.
- Estimating acceptable spiked species retrieval %
- Data transformation, Differential abundance and visualization.

*Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples based on our previous studies [WalkerLab](https://walkerlabmtsu.weebly.com/personnel.html).

---

## GCN Normalization with QIIME2 Plugin

Opinions on gene copy number (GCN) correction for the 16S rRNA marker vary, with proponents citing improved accuracy and critics noting limitations. While GCN correction is not included in the DspikeIn package, it can be applied to relative abundance counts using tools like the `q2-gcn-norm` plugin in Qiime2 (rrnDB v5.7) or methods outlined by [Louca et al., 2018](https://link.springer.com/content/pdf/10.1186/s40168-018-0420-9). Due to variability in rDNA copy numbers, GCN correction was not applied to ITS data.
We used the `q2-gcn-norm` plugin to normalize data by gene copy number (GCN). For more details, visit the [q2-gcn-norm GitHub repository](https://github.com/Jiung-Wen/q2-gcn-norm).

### Command Example

```bash
qiime gcn-norm copy-num-normalize \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --o-gcn-norm-table table-normalized.qza
```

---
*If you encounter issues installing the package due to missing dependencies, follow these steps to install all required packages first:*

## Step 1: Install Required Packages

To install the required packages, use the following script:

#### CRAN packages

```r
# Install CRAN packages
install.packages(c("stats", "dplyr", "ggplot2", "flextable","ggpubr", "randomForest", "ggridges", "ggalluvial","tibble", "matrixStats", "RColorBrewer", "ape", "rlang", "scales", "magrittr", "phangorn"))

# Load CRAN packages
lapply(c("stats", "dplyr", "ggplot2", "flextable","ggpubr","randomForest", "ggridges", "ggalluvial","tibble", "matrixStats", "RColorBrewer", "ape", "rlang", "scales", "magrittr", "phangorn"), library, character.only = TRUE)


```

#### Bioconductor Packages


```r
# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("phyloseq", "msa", "DESeq2","ggtree", "edgeR", "Biostrings", "DECIPHER", "microbiome"))

# Load Bioconductor packages
lapply(c("phyloseq", "msa", "DESeq2", "edgeR", "Biostrings","ggtree", "DECIPHER", "microbiome"), library, character.only = TRUE)

```

#### GitHub Packages


```r

# Install remotes if not installed
install.packages("remotes")

# Install GitHub packages
remotes::install_github("mikemc/speedyseq")
remotes::install_github("microsud/microbiomeutilities")

# Load GitHub packages
library(speedyseq)
library(microbiomeutilities)


```

## Step 2: Install DspikeIn Package



```r

# Installation
#Instructions for how to install the DspikeIn package.

# Using devtools
install.packages("devtools")
devtools::install_github("mghotbi/DspikeIn")
library(DspikeIn)

# Or using remotes
install.packages("remotes")
remotes::install_github("mghotbi/DspikeIn")
library(DspikeIn)

```


## Acknowledgement

DspikeIn builds on the excellent [**phyloseq**](https://github.com/joey711/phyloseq) package.

---



```r
# Make a new directory and set it as your working directory
create_directory("DspikeIn_16S_OTU", set_working_dir = TRUE)
getwd()


# Therefore, please start by creating a phyloseq object and follow the instructions.
# To create your phyloseq object, please refer to the phyloseq tutorial (https://joey711.github.io/phyloseq).
# The phyloseq object needs to include OTU/ASV, Taxa, phylogenetic tree, DNA reference, 
# and metadata containing spiked species volume, starting from 0 (no spike species added) to 4 (4 μl of spike cell added).

# Note: DspikeIn requires 'spiked.volume'; any other format is not readable."¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯"

# We are going to work with a subset of the dataset for both ASVs and OTUs
# approaches to accelerate this workshop.

physeq_16SOTU <-readRDS("physeq_16SOTU.rds")
physeq_ITSOTU <-readRDS("physeq_ITSOTU.rds")

physeq_16SOTU <- tidy_phyloseq(physeq_16SOTU)

# Ensure your metadata contains spiked volumes:
physeq_16SOTU@sam_data$spiked.volume


```

## Prepare the required information for our Protocol
#### Pre-process one Spiked-in Species


```r
# Required Information 
# Please note that the Spike cell numbers, species name, and selected hashcodes are customizable and can be tailored to the specific needs of individual studies.
# Moreover, to proceed with the DspikeIn package, you only need to select one method to specify your spiked species: either by hashcodes or species name.

library(phyloseq)
# 16S rRNA
# presence of 'spiked.volume' column in metadata
spiked_cells <-1847
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
merged_spiked_species<-"Tetragenococcus_halophilus"
Tetra <- subset_taxa(physeq_16SOTU,Species=="Tetragenococcus_halophilus" | Species=="Tetragenococcus_sp")
hashcodes <- row.names(phyloseq::tax_table(Tetra))

# ITS rDNA
# presence of 'spiked.volume' column in metadata
spiked_cells <- 733
species_name <- spiked_species<-merged_spiked_species<-"Dekkera_bruxellensis"
Dekkera <- subset_taxa(physeq_ITSOTU, Species=="Dekkera_bruxellensis")
hashcodes <- row.names(phyloseq::tax_table(Dekkera))

```

---

## Prepare the Required Information for the Synthetic Community
#### Pre-process a List of Spiked-in Species

---


```r

# Define the list of spiked-in species
spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")

# Define the corresponding copy numbers for each spiked-in species
spiked_cells_list <- c(10000, 20000, 15000) # or equal number of copies-> spiked_cells_list <- c(200, 200, 200)

```

## Plot phylogenetic tree with Bootstrap Values
This step will be helpful for handling ASVs with/without Gene Copy Number Correction
This section demonstrates how to use various functions from the package to plot and analyze phylogenetic trees.

---

```r
# In case there are several OTUs/ASVs resulting from the spiked species, you may want to check the phylogenetic distances.
# We first read DNA sequences from a FASTA file, to perform multiple sequence alignment and compute a distance matrix using the maximum likelihood method, then we construct a phylogenetic tree
# Use the Neighbor-Joining method  based on a Jukes-Cantor distance matrix and plot the tree with bootstrap values.
# we compare the Sanger read of Tetragenococcus halophilus with the FASTA sequence of Tetragenococcus halophilus from our phyloseq object.
# Load required libraries
  library(Biostrings)
  library(msa)
  library(phangorn)
  library(ape)
  library(speedyseq)
  library(ggtree)

# Subset the phyloseq object to include only Tetragenococcus species first
Tetra <- subset_taxa(Tetra, !is.na(taxa_names(Tetra)) & taxa_names(Tetra) != "")
tree <- phy_tree(Tetra)
ref_sequences_Tetra <- refseq(Tetra)
writeXStringSet(ref_sequences_Tetra, "ref_sequences_Tetra.fasta")
# postitive control 
Tetra_control_sequences <- Biostrings::readDNAStringSet("~/Tetra_Ju.fasta")

# combine the Tetragenococcus FASTA files (from your dataset and the Sanger fasta of Tetragenococcus, positive control)
combined_sequences <- c(ref_sequences_Tetra, Tetra_control_sequences)
writeXStringSet(combined_sequences, filepath = "~/combined_fasta_file")
combined_sequences <- Biostrings::readDNAStringSet("~/combined_fasta_file")

# Plot Neighbor-Joining tree with bootstrap values to compare Tetragenococcus in your dataset with your positive control
fasta_path <- "~/combined_fasta_file"
plot_tree_nj(fasta_path, output_file = "neighbor_joining_tree_with_bootstrap.png")


# Plot phylogenetic tree
plot_tree_custom(Tetra, output_prefix = "p0", width = 18, height = 18, layout = "circular")

# Plot the tree with glommed OTUs at 0.2 resolution/ or modify it
plot_glommed_tree(Tetra, resolution = 0.2, output_prefix = "top", width = 18, height = 18)

# Plot the phylogenetic tree with multiple sequence alignment
plot_tree_with_alignment(Tetra, output_prefix = "tree_alignment", width = 15, height = 15)

# Plot phylogenetic tree with bootstrap values and cophenetic distances
Bootstrap_phy_tree_with_cophenetic(Tetra, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 500)


```


Figure 1.

| Neighbor Joining Tree with Bootstrap | Tetra Plot with Bootstrap | Cophenetic tree with Bootstrap|
|:------------------------------------:|:----------:|:---------:|
| ![Neighbor Joining Tree with Bootstrap](https://github.com/mghotbi/DspikeIn/assets/29090547/175c8554-261f-45ab-8ec8-b22d09eb9ee4) | ![Tetra Plot with Bootstrap](https://github.com/mghotbi/DspikeIn/assets/29090547/78b5d8bb-3391-433d-af77-e40c8d27b055) | ![Cophenetic tree with Bootstrap](https://github.com/mghotbi/DspikeIn/assets/29090547/2a4bc232-e834-4212-9119-8561eddebed1) |




```markdown
## Aligned Sequences

The result of the aligned sequences is shown below:
DNAStringSet object of length 5:
    width seq                                                                                                                                         names               
[1]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGG-------------------- 2ddb215ff668b6a24...
[2]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGG-------------------- Tetragenococcus h...
[3]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTAGGTAGCAAACAGG-------------------- 65ab824f29da71010...
[4]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGACTCGAAAGCGTGGGTAGCAAACAGG-------------------- e49935179f23c00fb...
[5]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG-------------------- 0350f990080b4757a...

```


---

We selected the OTU approach using de novo robust clustering algorithms at a 97% similarity threshold, following the methods outlined by [Westcott and Schloss (2015)](https://doi.org/10.7717/peerj.1487).

Retrieved spiked species are correlated with abundance, richness, and beta diversity, and their recovery can therefore vary depending on the system under study. Our results indicate that the acceptable range of retrieved spiked species can be expanded to 35% in our model system. This contrasts with the findings of [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8), who reported an acceptable range of 0.1% to 10%.

---



| ASVs Or OTUs | Acceptable range|
|:----------:|:---------:|
| ![Desired range of spiked sp](https://github.com/mghotbi/DspikeIn/assets/29090547/2f949616-6493-4445-8f1e-7ac9c9dd844f) | ![Acceptable range](https://github.com/mghotbi/DspikeIn/assets/29090547/8674f3de-ba24-4857-9cd7-d1b6dc15c669) |



## Subsetting and Preprocessing Spiked Data


*Subset the part of the data which is spiked. Keep solely spiked samples using the `spiked.volume` column.*


```r

# Subset spiked samples (264 samples are spiked)
spiked_16S_OTU <- subset_samples(physeq16S_OTU, spiked.volume %in% c("2", "1"))
spiked_16S_OTU <- tidy_phyloseq(spiked_16S_OTU)

```

### Examine Your Count Table/Biom File Before Going Further


```r

# Summarize the initial statistics for ASVs/OTUs
initial_stat_ASV <- summ_phyloseq_ASV_OTUID(spiked_16S_OTU)

# Summarize the initial statistics sample-wise
initial_stat_sampleWise <- summ_phyloseq_sampleID(spiked_16S_OTU)

# Summarize the count data
summ_count_phyloseq(physeq_16S_OTU)

# Check the summary statistics
# Ensure the input is in dataframe format for this function
calculate_summary_stats_table(initial_stat_sampleWise)


```


*Check if transformation is required for spike volume variation.*


```r

# Adjust abundance by one-third
readAdj16S <- adjust_abundance_one_third(spiked_16S_OTU, factor = 3)
summ_count_phyloseq(readAdj16S)

# Random subsampling with reduction factor foe count and taxa
red16S <- random_subsample_WithReductionFactor(spiked_16S_OTU, reduction_factor = 3)
summ_count_phyloseq(red16S)


```


## Preprocessing for Scaling Factor Calculation  

If the spiked species appear in several OTUs/ASVs, check their phylogenetic distances and compare them to the reference sequences of your positive control.


```r

# Modify the threshold of acceptable spiked species % as needed. 
# For detailed guidance on acceptable thresholds (passed_range), 
# please refer to the instructions in our upcoming paper.

# Merge the spiked species
# merge_method = "max": Selects the maximum abundance among OTUs/ASVs of the spiked species, ensuring the most abundant ASV is retained.
# merge_method = "sum": Sums the abundances of OTUs/ASVs of the spiked species, providing a cumulative total.

species_name <- "Tetragenococcus_halophilus"

# Merge using "sum" or "max" method
Spiked_16S_sum_scaled <- Pre_processing_species(
  spiked_16S_OTU, 
  species_name, 
  merge_method = "sum", 
  output_file = "merged_physeq_sum.rds")

# Merge hashcodes using "sum" or "max" method
Spiked_16S_sum_scaled <- Pre_processing_hashcodes(
  spiked_16S_OTU, 
  hashcodes, 
  merge_method = "sum", 
  output_prefix = "merged_physeq_sum")

# Summarize count
summ_count_phyloseq(Spiked_16S_sum_scaled)

# Tidy phyloseq object
Spiked_16S_OTU_scaled <- tidy_phyloseq(Spiked_16S_sum_scaled)

# Now calculate the spiked species retrieval percentage.
# Customize the passed_range and merged_spiked_species/merged_spiked_hashcodes based on your preferences.
# passed_range = "c(0.1, 11)": threshold of acceptable spiked species %
# Select either merged_spiked_species or merged_spiked_hashcodes

merged_spiked_species <- c("Tetragenococcus_halophilus")
result <- calculate_spike_percentage(
  Spiked_16S_sum_scaled, 
  merged_spiked_species, 
  passed_range = c(0.1, 11))
calculate_summary_stats_table(result)

# Define your merged_spiked_hashcodes
merged_Tetra <- subset_taxa(
  Spiked_16S_OTU_scaled, 
  Species == "Tetragenococcus_halophilus")

merged_spiked_hashcodes <- row.names(tax_table(merged_Tetra))
result <- calculate_spike_percentage(
  Spiked_16S_OTU_scaled,  
  merged_spiked_hashcodes, 
  passed_range = c(0.1, 11))
calculate_summary_stats_table(result)

# If you decide to remove the failed reads and go forward with passed reads, here is what you need to do
# You can also go forward with the original file and remove the failed reads 
# after converting relative to absolute abundance

# Filter to get only the samples that passed
passed_samples <- result$Sample[result$Result == "passed"]

# Subset the original phyloseq object to keep only the samples that passed
passed_physeq <- prune_samples(
  passed_samples, 
  Spiked_16S_OTU_scaled)

```


### Estimating Scaling Factors After Pre-Processing

To estimate scaling factors, ensure you have the `merged_spiked_species` data, which contains the merged species derived from the spiking process.
*As we have already merged either hashcodes or spiked species and are aware of the contents of the taxa table, we can proceed from here with merged_spiked_species.*


```r
# Define the merged spiked species
merged_spiked_species <- c("Tetragenococcus_halophilus")

# Calculate spikeIn factors
result <- calculate_spikeIn_factors(Spiked_16S_OTU_scaled, spiked_cells, merged_spiked_species)

# Check the outputs
scaling_factors <- result$scaling_factors
physeq_no_spiked <- result$physeq_no_spiked
spiked_16S_total_reads <- result$spiked_16S_total_reads
spiked_species_reads <- result$spiked_species_reads

```


## Convert Relative Counts to Absolute Counts and Create a New Phyloseq Object


```r

# Convert relative counts data to absolute counts
absolute <- convert_to_absolute_counts(Spiked_16S_OTU_scaled, scaling_factors)
absolute_counts <- absolute$absolute_counts
physeq_absolute_abundance_16S_OTU <- absolute$physeq_obj


# summary statistics 
post_eval_summary <- calculate_summary_stats_table(absolute_counts)
print(post_eval_summary)


```


## Let's check the conclusion and get the report table of spiked species success or failure.



```r

# Define the parameters once.
merged_spiked_species <- c("Tetragenococcus_halophilus")
max_passed_range <- 35

# Subset the phyloseq object to exclude blanks
physeq_absolute_abundance_16S_OTU_perc <- subset_samples(physeq_absolute_abundance_16S_OTU, sample.or.blank != "blank")

# Generate the spike success report and summary statistics
summary_stats <- conclusion(physeq_absolute_abundance_16S_OTU_perc, merged_spiked_species, max_passed_range)
print(summary_stats)


```

Here is an example of a success or failure report:
![success report](https://github.com/mghotbi/DspikeIn/assets/29090547/017cfa65-8b75-4625-8d49-6e4a67146193)



```r

#Save your file for later. Please stay tuned for the rest: Comparisons and several visualization methods to show how important it is to convert relative to absolute abundance in the context of microbial ecology.

physeq_absolute_16S_OTU <- tidy_phyloseq(physeq_absolute_abundance_16S_OTU_perc)
saveRDS(physeq_absolute_16S_OTU, "physeq_absolute_16S_OTU.rds")

```
## Normalization and bias correction 


```r
# Bolstad, B.M., Irizarry, R.A., Åstrand, M. and Speed, T.P., 2003. A comparison of normalization methods for high-density oligonucleotide array data based on variance and bias. Bioinformatics, 19(2), pp.185-193.
# Gagnon-Bartsch, J.A. and Speed, T.P., 2012. Using control genes to correct for unwanted variation in microarray data. Biostatistics, 13(3), pp.539-552.
# Risso, D., Ngai, J., Speed, T.P. and Dudoit, S., 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature biotechnology, 32(9), pp.896-902.
# Gagnon-Bartsch, J.A., Jacob, L. and Speed, T.P., 2013. Removing unwanted variation from high dimensional data with negative controls. Berkeley: Tech Reports from Dep Stat Univ California, pp.1-112.

# Load required libraries
library(phyloseq)
library(DESeq2)
library(edgeR)
library(BiocGenerics)

#ps is a phyloseq object without spiked species counts
ps <- physeq_absolute_16S_OTU
ps <- remove_zero_negative_count_samples(physeq_absolute_abundance_16S_OTU)
ps <- convert_categorical_to_factors(physeq_absolute_abundance_16S_OTU)
# group_var <- "Animal.ecomode"  

# Normalization Methods:
# result_DESeq <- normalization_set(ps, method = "DESeq", groups = "group_var")
# result_TMM <- normalization_set(ps, method = "TMM", groups = "group_var")
# result_CLR <- normalization_set(ps, method = "clr")
# result_UQ <- normalization_set(ps, method = "UQ", groups = group_var)
# result_med <- normalization_set(ps, method = "med", groups = group_var)
# result_Poisson <- normalization_set(ps, method = "Poisson", groups = "group_var")
# result_UQ <- normalization_set(ps, method = "UQ", groups = "group_var")
# result_med <- normalization_set(ps, method = "med", groups = "group_var")
# result_rle <- normalization_set(ps, method = "rle")
# result_css <- normalization_set(ps, method = "CSS")
# result_tss <- normalization_set(ps, method = "tss")
# result_rar <- normalization_set(ps, method = "rar")

# Customized filtering and transformations
# Proportion adjustment
normalized_physeq <- proportion_adj(ps, output_file = "proportion_adjusted_physeq.rds")
summ_count_phyloseq(normalized_16S)


# Relativize and filter taxa based on selected thresholds
FT_physeq <- relativized_filtered_taxa(
  ps,
  threshold_percentage = 0.0001,
  threshold_mean_abundance = 0.0001,
  threshold_count = 5,
  threshold_relative_abundance = 0.0001)
summ_count_phyloseq(FT_physeq)

# Adjust prevalence based on the minimum reads
physeq_min <- adjusted_prevalence(ps, method = "min")


```


---

## Visualization and Differential abundance 


```r


# taxa barplot
#abundance_type = "absolute"/"relative"
bp_ab <- taxa_barplot(physeq_absolute_16S_OTU, target_glom = "Genus", treatment_variable = "Host.genus", abundance_type = "absolute", x_angle = 90, fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
print(bp_ab$barplot)

# original relative count -> spiked_16S_OTU
bp_rel <- taxa_barplot(spiked_16S_OTU, target_glom = "Genus", treatment_variable = "Host.genus", abundance_type = "relative", x_angle = 90, fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
print(bp_rel$barplot)

```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Rel AbsSal taxa barplot](https://github.com/mghotbi/DspikeIn/assets/29090547/643fca2a-6087-49f1-b3ee-5759d2fcb36f) | ![Rel abun Sal taxa barplot](https://github.com/mghotbi/DspikeIn/assets/29090547/2040830e-1ce1-46c9-8dfe-3c18f77a85bf) |




```r

# simple barplot of taxonomy abundance
# Plot relativized abundance
plot <- plotbar_abundance(physeq_absolute_16S_OTU, level = "Family", group = "Env.broad.scale.x", top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = TRUE, output_prefix = "relativized_abundance_plot")
print(plot)

# Plot non-relativized (absolute) abundance
plot_absolute <- plotbar_abundance(spiked_16S_OTU, level = "Family", group = "Env.broad.scale.x", top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE, output_prefix = "non_relativized_abundance_plot")
print(plot_absolute)


# Check abundance distribution via Ridge Plots before and after converting to absolute abundance
ridgeP_before <- ridge_plot_it(spiked_16S_OTU, taxrank = "Family", top_n = 10)
ridgeP_after <- ridge_plot_it(physeq_absolute_16S_OTU, taxrank = "Family", top_n = 10)


```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/7b50556d-77d7-4aa0-b2b0-f461c67c65a7) | ![Rel  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/5b5a85ea-1f10-4082-a108-7c9019bd84d8) |


```r

# core_microbiome
custom_detections <- 10^seq(log10(3e-1), log10(0.5), length = 5)
PCM_rel <- plot_core_microbiome_custom(spiked_16S_OTU, detections = custom_detections, taxrank = "Family", output_core_rds = "core_microbiome.rds", output_core_csv = "core_microbiome.csv")

PCM_Abs <- plot_core_microbiome_custom(physeq_absolute_16S_OTU, detections = custom_detections, taxrank = "Family", output_core_rds = "core_microbiome.rds", output_core_csv = "core_microbiome.csv")

# core.microbiome is automatically saved in your working directory so yoou can go ahead and barplot it
core.microbiome <- readRDS("core.microbiome.rds")

```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs core](https://github.com/mghotbi/DspikeIn/assets/29090547/9f32f799-4421-4842-a4cf-a1eee1a768e1) | ![Rel core](https://github.com/mghotbi/DspikeIn/assets/29090547/5fb4055a-33f5-48ec-bb0c-a9a9ab446429) |


```r

# shift to long-format data frame and plot the abundance of taxa across the factor of your interest
# Generate alluvial plot
ps_physeq_absolute_16S_OTU <- psmelt(physeq_absolute_16S_OTU)

# Define total reads for relative abundance calculation
total_reads <- sum(ps_physeq_absolute_16S_OTU$Abundance)  
# Generate alluvial plot for absolute abundance
alluvial_plot_abs <- alluvial_plot(data = ps_physeq_absolute_16S_OTU,
                                   axes = c("Host.genus", "Ecoregion.III"),
                                   abundance_threshold = 1000, fill_variable = "Family",
                                   silent = TRUE, abundance_type = "absolute",
                                   top_taxa = 10, facet_vars = c("Diet"))



```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs Alluv](https://github.com/mghotbi/DspikeIn/assets/29090547/2f187727-db7b-41a2-82be-73162423ce25) | ![Rel Alluv](https://github.com/mghotbi/DspikeIn/assets/29090547/bc6ed255-97d3-4e24-ad22-12890b747e79) |


```r

# selecting the most important ASVs/OTUs through RandomForest classification
# Salamander_absolute= subset of our phyloseq object
rf_physeq <- RandomForest_selected_ASVs(ps_physeq_absolute_16S_OTU, response_var = "Host_Species", na_vars = c("Habitat","Diet", "Ecoregion_III", "Host_genus", "Animal_type"))
RP=ridge_plot_it(rf_physeq)
RP+facet_wrap(~Diet)


#detect common ASVs/OTUs
# The input is the list of phyloseq objects
results <- detect_common_asvs_taxa(list(rf_physeq, FTspiked_16S , core.microbiome), 
                                    output_common_asvs_rds = "common_asvs.rds", 
                                    output_common_taxa_rds = "common_taxa.rds")

common_asvs_phyloseq <- results$common_asvs_phyloseq
common_taxa_phyloseq <- results$common_taxa_phyloseq

plotbar_abundance(common_taxa_phyloseq, level = "Family", group = "Env.broad.scale", top = 10, return = TRUE)

```
