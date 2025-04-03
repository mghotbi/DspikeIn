## Welcome to the DspikeIn R package repository!

**Author:** Mitra Ghotbi  
 **Version:** 1.2.3  
 **Date:** March 8, 2025  



![CheatSheetDspikeIn](https://github.com/user-attachments/assets/1972dbe7-05b6-4d61-a683-d0ac02b2a219)

---

## 📚 Table of Contents

1. **Getting Started**
   DspikeIn Package
   - [DspikeIn](#DspikeIn-Package)
   - [Installation](#Installation)
   - [Requirements](#To-Meet-Taxonomic-Ranks-Requirements)
   - [GCN Correction](#GCN-Normalization-with-QIIME2-Plugin)
   - [Dataset for training](#dataset-for-practicing-dspikein-package)
   - [Build Phyloseq or TSE](#Building-your-own-phyloseq-and-TSE)



3. **Data Preparation**
   - [Validation Using Phylogenetic Tree](#spike-in-validation)
   - [Spike-ins Behaviour ](#Spike-ins-Behaviour)
   - [Preparation for Our Protocol](#Prepare-the-required-information-for-our-protocol)
   - [Preparation for the Synthetic Community](#Prepare-the-Required-Information-for-the-Synthetic-community)

4. **Processing**
   - [Preprocessing One Species Scaling Factor](#Preprocessing-One-Species-Scaling-Factor)
   - [Preprocessing List of Species Scaling Factor](#Preprocessing-List-of-Species-Scaling-Factor)
   - [Calculate Spiked Species Retrieval % for One Species](#calculate-spiked-species-retrieval--for-one-species)
   - [Calculate Spiked Species Retrieval % for List of Species](#calculate-spiked-species-retrieval--for-list-of-species)
   - [Scaling Factors for One Spiked Species](#scaling-factors-for-one-spiked-species)
   - [Scaling Factors for a List of Spiked Species](#scaling-factors-for-a-list-of-spiked-species)
   - [System-Specific Spiked Species Retrieval](#System-specific-spiked-species-retrieval)
   - [Conclusion](#conclusion)


5. **Bias Correction**
   - [Convert Relative to Absolute Counts](#convert-relative-counts-to-absolute-counts-and-create-a-new-phyloseq-object)
   - [Normalization and Differential Abundance](#normalization-and-differential-abundance)
   - [Customized Filtering](#customized-filtering)

6. **Visualization**
   - [Visualization](#visualization)
   - [Detect common ASVs/OTUs](#Detect-common-asvs-otus)

7. **Credits**
   - [Acknowledgement](#Acknowledgement)
   - [Citing DspikeIn](#if-you-use-this-package-and-find-it-useful)

---

### DspikeIn Package

DspikeIn is designed for microbiome data analysis, seamlessly integrating with **phyloseq** (for marker-gene microbiome data) and **TreeSummarizedExperiment (TSE)** (for hierarchical biological data, including microbiomes). These objects must include seven taxonomic ranks.  
For absolute abundance estimation, the metadata must contain **spiked.volume**.  

## Features of DspikeIn  
The DspikeIn package provides functions for:  

- Verifying the phylogenetic distances of ASVs/OTUs derived from spiked species.  
- Preprocessing microbiome data.  
- Calculating spike-in scaling factors.  
- Converting relative abundance to absolute abundance.  
- Estimating acceptable retrieval percentages of spiked species.  
- Performing data transformation, differential abundance analysis, and visualization.  


### To get detailed examples and guidance, please use:
browseVignettes("DspikeIn")

### Data availability
The DspikeIn package provides example datasets located in the data/ folder and inst/extdata/ folder. You can list the available datasets using the following commands:

```r

# List datasets available in the DspikeIn package
data(package = "DspikeIn")

# List files in the extdata folder
list.files(system.file("extdata", package = "DspikeIn"))
```
### Building your own phyloseq and TSE

```r

# =====================================================================
#                     Build phyloseq 
# =====================================================================
otu <- read.csv("otu.csv", header = TRUE, sep = ",", row.names = 1)
# taxonomic rank need to be capilalized, only the first letter of each rank
tax <- read.csv("tax.csv", header = TRUE, sep = ",", row.names = 1)
# Ensure 'spiked.volume' column is present and correctly formatted in metadata
meta <- read.csv("metadata.csv", header = TRUE, sep = ",")

# Convert data to appropriate formats
meta <- as.data.frame(meta)
taxmat <- as.matrix(tax)
otumat <- as.matrix(otu)
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxmat)

# Check
row.names(meta) <- sample_names(OTU)
metadata <- sample_data(meta)
# Build phyloseq obj
physeq <- phyloseq(OTU, TAX, metadata)

# Follow the next steps if tree and reference files are included
MyTree <- read.tree("tree.nwk")
reference_seqs <- readDNAStringSet(file = "dna-sequences.fasta", format = "fasta")

physeq_16SOTU <- merge_phyloseq(physeq, reference_seqs, MyTree)
physeq_16SOTU <- tidy_phyloseq_tse(physeq_16SOTU)

saveRDS(physeq_16SOTU, file = "physeq_16SOTU.rds")
physeq_16SOTU <- readRDS("physeq_16SOTU.rds")


# =====================================================================
#                       Build TSE 
# =====================================================================

otu <- read.csv("otu.csv", header = TRUE, sep = ",", row.names = 1)
otu_mat <- as.matrix(otu)  # Convert to matrix
tax <- read.csv("tax.csv", header = TRUE, sep = ",", row.names = 1)
colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")  
tax_mat <- as.matrix(tax)  # Convert to matrix
meta <- read.csv("metadata.csv", header = TRUE, sep = ",", row.names = 1)
reference_seqs <- readDNAStringSet("dna-sequences.fasta", format = "fasta")
tse <- TreeSummarizedExperiment(
  assays = list(counts = otu_mat),  # OTU table 
  rowData = tax_mat,                # Taxonomy information
  colData = meta,                    # Sample metadata
  rowTree = MyTree,                  # Phylogenetic tree
  rowSeqs = reference_seqs           # Reference sequences
)


```
**Whole-Cell Spike-In Protocol,**
*Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples based on our previous studies [WalkerLab](https://walkerlabmtsu.weebly.com/personnel.html).

---
### GCN Normalization with QIIME2 Plugin

Opinions on gene copy number (GCN) correction for the 16S rRNA marker vary, with proponents citing improved accuracy and critics noting limitations. While GCN correction is not included in the DspikeIn package, it can be applied to relative abundance counts using tools like the `q2-gcn-norm` plugin in Qiime2 (rrnDB v5.7) or methods outlined by [Louca et al., 2018](https://link.springer.com/content/pdf/10.1186/s40168-018-0420-9),including PICRUSt, CopyRighter, and PAPRICA. Due to variability in rDNA gene copy numbers ([Lavrinienko et al., 2021](https://doi.org/10.1186/s42523-021-00134-z)), GCN corrections were not applied. However, targeted adjustments can be made to prevent overestimating specific fungal taxa.

---
### Command Example

```bash
qiime gcn-norm copy-num-normalize \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --o-gcn-norm-table table-normalized.qza

```

## Installation

To install the DspikeIn package, follow these steps...
*If you encounter issues installing the package due to missing dependencies, follow these steps to install all required packages first:*

## Step 1: Install Required Packages
To install the required packages, use the following script:

---
#### CRAN packages

```r
# Install missing CRAN packages
install.packages(setdiff(c("stats", "dplyr", "ggplot2", "flextable", "ggpubr", 
                           "randomForest", "ggridges", "ggalluvial", "tibble", 
                           "matrixStats", "RColorBrewer", "ape", "rlang", 
                           "scales", "magrittr", "phangorn", "igraph", "tidyr", 
                           "xml2", "data.table", "reshape2","vegan", "patchwork", "officer"), 
                         installed.packages()[,"Package"]))

# Load CRAN packages
lapply(c("stats", "dplyr", "ggplot2", "flextable", "ggpubr", "randomForest", 
         "ggridges", "ggalluvial", "tibble", "matrixStats", "RColorBrewer", 
         "ape", "rlang", "scales", "magrittr", "phangorn", "igraph", "tidyr", 
         "xml2", "data.table", "reshape2","vegan", "patchwork", "officer"), library, character.only = TRUE)

```
#### Bioconductor Packages

```r 

# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install missing Bioconductor packages
BiocManager::install(setdiff(c("phyloseq", "msa", "DESeq2", "ggtree", "edgeR", 
                               "Biostrings", "DECIPHER", "microbiome", "limma", 
                               "S4Vectors", "SummarizedExperiment", "TreeSummarizedExperiment"), 
                             installed.packages()[,"Package"]))

# Load Bioconductor packages
lapply(c("phyloseq", "msa", "DESeq2", "edgeR", "Biostrings", "ggtree", "DECIPHER", 
         "microbiome", "limma", "S4Vectors", "SummarizedExperiment", "TreeSummarizedExperiment"), 
       library, character.only = TRUE)

```
#### GitHub Packages

```r

# Install remotes if not installed
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
library(remotes)

# Install missing GitHub packages
remotes::install_github("mikemc/speedyseq")
remotes::install_github("microsud/microbiomeutilities")
# Optional
#devtools::install_github("briatte/ggnet")
#devtools::install_github("zdk123/SpiecEasi")

# Load GitHub packages
library(speedyseq)
library(microbiomeutilities)
#library(SpiecEasi)
#library(ggnet)


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

# To access the DspikeIn vignette for a detailed tutorial, use vignette("DspikeIn"), or browse all available vignettes with browseVignettes("DspikeIn").
devtools::install_github("mghotbi/DspikeIn", build_vignettes = TRUE, dependencies = TRUE)
browseVignettes("DspikeIn")
vignette("DspikeIn")

```


## Acknowledgement

DspikeIn builds on the excellent [**phyloseq**](https://github.com/joey711/phyloseq) package.
Requirements


### To Meet Taxonomic Ranks Requirements

```r

### To remove strain from the taxonomic ranks

library(phyloseq)
# Function to remove strain information from taxonomy columns
remove_strain_info <- function(tax_table) {
  # Define the regex pattern for common strain identifiers
  pattern <- "Strain.*|strain.*|\\s*\\[.*\\]|\\s*\\(.*\\)"  # Adjust the regex to match specific strain formats
  # Apply the pattern to each column of the taxonomy table
  for (col in colnames(tax_table)) { 
    tax_table[, col] <- gsub(pattern, "", tax_table[, col])  # Remove strain info
    tax_table[, col] <- trimws(tax_table[, col])  # Trim trailing whitespace
  }
    return(tax_table)}
# Step 1: Extract the taxonomy table
taxonomy <- tax_table(ps)
# Step 2: Remove strain information (including the `Strain` column)
cleaned_taxonomy <- remove_strain_info(taxonomy)
# Remove the `Strain` column if it exists
if ("Strain" %in% colnames(cleaned_taxonomy)) {
  cleaned_taxonomy <- cleaned_taxonomy[, colnames(cleaned_taxonomy) != "Strain"]}
# Step 3: Update the taxonomy table in the phyloseq object
tax_table(ps) <- cleaned_taxonomy
# Step 4: Verify the changes
print(head(tax_table(ps)))  # Display the first few rows

```

To add species rank to the taxonomic ranks


```r
library(phyloseq)
# Step 1: Extract the taxonomy table from the phyloseq object
taxonomy <- tax_table(ps)
# Step 2: Add a `species` column based on the Genus column and unique identifiers
taxonomy[, "species"] <- paste0(taxonomy[, "Genus"], "_OTU", seq_len(nrow(taxonomy)))
# Step 3: Update the taxonomy table in the phyloseq object
tax_table(ps) <- taxonomy
# Step 4: Verify the changes
print(head(tax_table(ps)))  


```

## Dataset for practicing DspikeIn Package

You can download the practice dataset for the DspikeIn package by clicking the link below:  
👉 [Download Dataset](https://drive.google.com/drive/folders/164_K7MaFLCf5T8F9fsPAb1AndQ8mJOOP?usp=sharing)

### OR You can directly use the datasets provided in the DspikeIn package to practice and test the functions. These datasets include phyloseq objects for different marker-gene microbiome analyses:

```r
# Load 16S rRNA OTU dataset
data("physeq_16SOTU", package = "DspikeIn")  

# Load ITS OTU dataset (for fungal microbiome analysis)
data("physeq_ITSOTU", package = "DspikeIn")

```

#### Make a new directory and set it as your working directory
```r

create_directory("DspikeIn_16S_OTU", set_working_dir = TRUE)
getwd()

# Therefore, please start by creating a phyloseq object and follow the instructions.
# To create your phyloseq object, please refer to the phyloseq tutorial (https://joey711.github.io/phyloseq).
# The phyloseq object needs to include OTU/ASV, Taxa, phylogenetic tree, DNA reference, 
# and metadata containing spiked species volume, starting from 0 (no spike species added) to 4 (4 μl of spike cell added).

# Note: DspikeIn requires 'spiked.volume'; any other format is not readable."¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯"

# We will work with a subset of the dataset for both ASV and OTU approaches to accelerate this workshop. However, the full dataset with the OTU approach is available
# in DspikeIn via data("physeq_16SOTU", package = "DspikeIn") and data("physeq_ITSOTU", package = "DspikeIn").

physeq_16SOTU <-readRDS("Relative16SOTU.rds")
physeq_ITSOTU <-readRDS("RelativeITSOTU.rds")

physeq_16SOTU <- tidy_phyloseq_tse(physeq_16SOTU)

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

## Prepare the Required Information for the Synthetic Community
 Pre-process List of Spiked-in Species

```r

# Define the list of spiked-in species
spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")

# Define the corresponding copy numbers for each spiked-in species
spiked_cells_list <- c(10000, 20000, 15000) # or equal number of copies-> spiked_cells_list <- c(200, 200, 200)

```

Do all detected sample spike-in sequences cluster with the reference, and are their branch lengths statistically similar, supporting a common ancestor?

# spike-in validation 

All sample-derived sequences are forming a clade with the reference.
We look for a monophyletic grouping of spike-in OTUs
The clade is strongly supported (bootstrap around 100 percentage).
The branch lengths and distances are in a biologically plausible range.

```r

# Use the Neighbor-Joining method  based on a Jukes-Cantor distance matrix and plot the tree with bootstrap values.
# we compare the Sanger read of Tetragenococcus halophilus with the FASTA sequence of Tetragenococcus halophilus from our phyloseq object.

library(Biostrings)
library(phyloseq)

# Get path to external data folder
extdata_path <- system.file("extdata", package = "DspikeIn")
list.files(extdata_path)

data("physeq_16SOTU", package = "DspikeIn")

physeq_16SOTU <- tidy_phyloseq_tse(physeq_16SOTU)

physeq_16SOTU <- phyloseq::prune_taxa(
  get_tax_table(physeq_16SOTU)$Kingdom %in% c("Bacteria", "Archaea"),
  physeq_16SOTU
)


# Subset the phyloseq object to include only Tetragenococcus species first
# Tetragenococcus <- subset_taxa(physeq_16SOTU, Genus=="Tetragenococcus")
# Tetragenococcus <- subset_taxa(Tetragenococcus, !is.na(taxa_names(Tetragenococcus)) & taxa_names(Tetragenococcus) != "")
# tree <- phy_tree(Tetragenococcus)
# ref_sequences_Tetragenococcus <- refseq(Tetragenococcus)
# library(Biostrings)
# writeXStringSet(ref_sequences_Tetragenococcus, "Sample.fasta.fasta")



ref_fasta <- system.file("extdata", "Ref.fasta", package = "DspikeIn")
sample_fasta <- system.file("extdata", "Sample.fasta", package = "DspikeIn")



result <- validate_spikein_clade(
  reference_fasta = ref_fasta,
  sample_fasta = sample_fasta,
  bootstrap = 200,
  output_prefix = "comp_spikein")

```

| Spike-in distribution | Spike-in validation |
|:---------------------:|:-------------------:|
| <img src="https://github.com/user-attachments/assets/a883a66a-9155-4a92-be5a-996cc4f2315e" width="500"/> | <img src="https://github.com/user-attachments/assets/406cacf3-fb06-480a-8a35-1cd800a92b57" width="500"/> |


### Spike-ins Behaviour 
####     Did spike-ins behave as expected across all samples?

```r

Tip labels=	OTU/ASV names
Branch length numbers=	Actual evolutionary distances (small = very similar)
Prevalence stars	How frequently the OTU occurs across samples
Blue bar ring=	Log10 mean abundance
Outer colored tiles=	The metadata variable you choose (e.g., Animal.type)

spikein <- phyloseq::subset_taxa(physeq_16SOTU, Genus == "Tetragenococcus")
taxa_names(spikein) <- paste0("OTU", seq_len(ntaxa(spikein)))

# Visualize
ps <- plot_spikein_tree_diagnostic(
  obj = spikein,
  metadata_var = "Animal.type",
  output_prefix = "tetragenococcus_diag")

ps
                                                                                                                               
```

We selected the OTU approach using de novo robust clustering algorithms at a 97% similarity threshold, following the methods outlined by [Westcott and Schloss (2015)](https://doi.org/10.7717/peerj.1487).

## Subsetting and Preprocessing Spiked Data
*Subset the part of the data which is spiked. Keep solely spiked samples using the `spiked.volume` column.*


```r
# Subset spiked samples (264 samples are spiked)
spiked_16S_OTU <- subset_samples(physeq16S_OTU, spiked.volume %in% c("2", "1")) # Optional
spiked_16S_OTU <- tidy_phyloseq_tse(spiked_16S_OTU)

```

Examine Your Count Table/Biom File Before Going Further

```r

# Summarize the initial statistics sample-wise
initial_stat_sampleWise <- summ_phyloseq_sampleID(spiked_16S_OTU)

# Summarize the count data
summ_count_phyloseq(spiked_16S_OTU)

# Check the summary statistics
# Ensure the input is in dataframe format for this function
calculate_summary_stats_table(initial_stat_sampleWise)

```

*Check if transformation is required for spike volume variation.*

```r

# Random subsampling with reduction factor foe count and taxa
red16S <- random_subsample_WithReductionFactor(spiked_16S_OTU, reduction_factor = 3)
summ_count_phyloseq(red16S)


```

## Preprocessing before Scaling Factor Calculation  
### Preprocessing One Species Scaling Factor
 
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

# Merge hashcodes using "sum" or "max" method ##for QIIME users
Spiked_16S_sum_scaled <- Pre_processing_hashcodes(
  spiked_16S_OTU, 
  hashcodes, 
  merge_method = "sum", 
  output_prefix = "merged_physeq_sum")

# Summarize count
summ_count_phyloseq(Spiked_16S_sum_scaled)

# Tidy phyloseq object
Spiked_16S_OTU_scaled <- tidy_phyloseq_tse(Spiked_16S_sum_scaled)


# Customize the passed_range and merged_spiked_species/merged_spiked_hashcodes based on your preferences.

```

### Spiked Species Retrieval % 
### Calculate Spiked Species Retrieval % for One Species

```r

# Customize the passed_range and merged_spiked_species/merged_spiked_hashcodes based on your preferences.
# passed_range = "c(0.1, 11)": threshold of acceptable spiked species %
# Select either merged_spiked_species or merged_spiked_hashcodes

merged_spiked_species <- c("Tetragenococcus_halophilus")

result <- calculate_spike_percentage(
  Spiked_16S_sum_scaled, 
  merged_spiked_species,
  output_file = "merged_perc_data.docx", 
  passed_range = c(0.1, 11))


calculate_summary_stats_table(result)

# result of calculating spiked sp reterival can be read as csv
merged_perc_data <- read.csv("merged_perc_data.csv")


## For QIIME users
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

# Two next steps are optional
# Filter to get only the samples that passed
passed_samples <- result$Sample[result$Result == "passed"]

# Subset the original phyloseq object to keep only the samples that passed
passed_physeq <- prune_samples(
  passed_samples, 
  Spiked_16S_OTU_scaled)

```

### Calculate Spiked Species Retrieval % for List of Species

```r

#ps= is phyloseq obj
spiked_species_list <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
result <- calculate_spike_percentage_list(ps, merged_spiked_species = spiked_species_list, passed_range = c(0.1, 10)) #change the range based on your system specifications
print(result)

```

### Preprocessing List of Species Scaling Factor  

```r

spiked_species <- c("Pseudomonas aeruginosa", "Escherichia coli", "Clostridium difficile")
merged_physeq_sum <- Pre_processing_species_list(physeq, spiked_species, merge_method = "sum")


```

### system specific spiked species retrieval
### Estimating the system-specific optimal range of spiked species retrieval using biological metrics
Previously, [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8) reported an acceptable range of 0.1% to 10% for spiked species retrieval. Here, we demonstrate that retrieved spiked species are correlated with biological metrics such as abundance, richness, evenness, and beta dispersion, indicating that the acceptable range is system-dependent and can be changed based on system specifications.


```r

# Results of Spiked Species Retrieval Can Be Added to Metadata
merged_perc_data <- read.csv("merged_perc_data.csv")

filtered_sample_data <- microbiome::meta(physeq_absolute) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Sample") %>%  
  dplyr::mutate(Sample = as.character(Sample)) %>%
  dplyr::left_join(merged_perc_data, by = "Sample")

filtered_sample_data <- tibble::column_to_rownames(filtered_sample_data, "Sample")

filtered_sample_data <- sample_data(as.data.frame(filtered_sample_data))

# Assign back to phyloseq obj
sample_data(physeq_absolute) <- filtered_sample_data

library(vegan)

# Calculate Pielou's Evenness using Shannon index and species richness (Observed)
alphab <- estimate_richness(physeq_absolute, measures = c("Observed", "Shannon"))
alphab$Pielou_evenness <- alphab$Shannon / alphab$Observed

# Normalize values
alphab <- alphab %>%
  mutate(across(c("Observed", "Shannon", "Pielou_evenness"), ~ as.numeric(scale(.))))

metadata <- as.data.frame(microbiome::meta(physeq_absolute))

metadata$Sample <- rownames(metadata)
alphab$Sample <- rownames(alphab)

# Merge alpha diversity metrics into metadata
metadata <- dplyr::left_join(metadata, alphab[, c("Sample", "Observed", "Shannon", "Pielou_evenness")], by = "Sample")

metadata <- metadata %>%
  column_to_rownames(var = "Sample")

# Updated metadata back to the phyloseq obj
sample_data(physeq_absolute) <- sample_data(metadata)

if (!"Spiked_Reads" %in% colnames(metadata)) {
  stop("Column 'Spiked_Reads' not found in metadata.")
}

# Generate regression plot
plot_object <- regression_plot(
  data = metadata,
  x_var = "Pielou_evenness",
  y_var = "Spiked_Reads",
  custom_range = c(0.1, 15, 30, 45, 60, 85, 100),
  plot_title = NULL
)

```

| Beta Dispersion  (16S) | Evenness (16S) |
|:---------------------:|:--------------:|
| ![BetaDis16S](https://github.com/user-attachments/assets/23178086-e7d3-4b0c-873d-5aefe0a12a8d) | ![Evenness16S](https://github.com/user-attachments/assets/3fc63f6f-50ea-4832-983f-68cd8cdf25a8) |

---

## Estimating Scaling Factors After Pre-Processing

### Scaling Factors for One Spiked Species

To estimate scaling factors, ensure you have the `merged_spiked_species` data, which contains the merged species derived from the spiking process.
*As we have already merged either hashcodes or spiked species and are aware of the contents of the taxa table, we can proceed from here with merged_spiked_species.*

```r
# Define the merged spiked species
merged_spiked_species <- c("Tetragenococcus_halophilus")

# Calculate spikeIn factors
result <- calculate_spikeIn_factors(Spiked_16S_OTU_scaled, spiked_cells, merged_spiked_species)

# Check the outputs
scaling_factors <- result$scaling_factors

```

### Scaling Factors for a List of Spiked Species

This example demonstrates how to calculate scaling factors after merging redundant spike-in species

```r
# Define spiked-in species and their corresponding cell counts
spiked_species_list <- list(
  c("Pseudomonas aeruginosa"),
  c("Escherichia coli"),
  c("Clostridium difficile"))

spiked_cells_list <- c(10000, 20000, 15000)

scaling_factors <- calculate_list_average_scaling_factors(
  merged_physeq_sum,
  spiked_species_list,
  spiked_cells_list,
  merge_method = "sum" ) # or max

```

## Convert Relative Counts to Absolute Counts and Create a New Phyloseq Object

```r
# Convert relative counts data to absolute counts
absolute <- convert_to_absolute_counts(Spiked_16S_OTU_scaled, scaling_factors)
absolute_counts <- absolute$absolute_counts
physeq_absolute_abundance_16S_OTU <- absolute$obj_adj

physeq_absolute <- physeq_absolute_abundance_16S_OTU

# summary statistics 
post_eval_summary <- calculate_summary_stats_table(absolute_counts)
print(post_eval_summary)


```

## Let's check the conclusion and get the report table of spiked species success or failure.
### Conclusion

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

physeq_absolute_16S_OTU <- tidy_phyloseq_tse(physeq_absolute_abundance_16S_OTU_perc)
saveRDS(physeq_absolute_16S_OTU, "physeq_absolute_16S_OTU.rds")

```

## Normalization and Differential Abundance

```r
# Bolstad, B.M., Irizarry, R.A., Åstrand, M. and Speed, T.P., 2003. A comparison of normalization methods for high-density oligonucleotide array data based on variance and bias. Bioinformatics, 19(2), pp.185-193.
# Gagnon-Bartsch, J.A. and Speed, T.P., 2012. Using control genes to correct for unwanted variation in microarray data. Biostatistics, 13(3), pp.539-552.
# Risso, D., Ngai, J., Speed, T.P. and Dudoit, S., 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature biotechnology, 32(9), pp.896-902.
# Gagnon-Bartsch, J.A., Jacob, L. and Speed, T.P., 2013. Removing unwanted variation from high dimensional data with negative controls. Berkeley: Tech Reports from Dep Stat Univ California, pp.1-112.


# Normalization Methods:
# group_var <- "Diet"  
# result_TMM <- normalization_set(ps, method = "TMM", groups = "group_var")
# result_UQ <- normalization_set(ps, method = "UQ", groups = group_var)
# result_med <- normalization_set(ps, method = "med", groups = group_var)
# result_Pois <- normalization_set(ps, method = "Poisson", groups = "group_var")
# result_DES <- normalization_set(ps, method = "DESeq", groups = "group_var")
# result_med <- normalization_set(ps, method = "med", groups = "group_var")
# result_rle <- normalization_set(ps, method = "rle")
# result_css <- normalization_set(ps, method = "CSS")
# result_tss <- normalization_set(ps, method = "tss")
# result_rar <- normalization_set(ps, method = "rar")
# result_CLR <- normalization_set(ps, method = "clr")

# lets subset samples we need from absolute data
# Absolute count
absolute_Des <- physeq_absolute %>%
  phyloseq::subset_taxa(Genus != "Tetragenococcus") %>%
  phyloseq::subset_samples(Clade.Order == "Caudate") %>%
  phyloseq::subset_samples(Host.genus %in% c("Desmognathus", "Plethodon", "Eurycea")) %>%
  phyloseq::subset_samples(Ecoregion.III == "Blue Ridge") %>%
  phyloseq::subset_samples(Host.genus == "Desmognathus")

results_DESeq2 <- perform_and_visualize_DA(
  obj = absolute_Des,
  method = "DESeq2",
  group_var = "Host.taxon",
  contrast = c("Desmognathus monticola", "Desmognathus imitator"),
  output_csv_path = "DA_DESeq2.csv",
  target_glom = "Genus",
  significance_level = 0.05
)
# Visualization
print(results_DESeq2$plot)
# Significant taxa table
head(results_DESeq2$results)
# Filtered phyloseq object with significant taxa
results_DESeq2$obj_significant


# Relative 
data("physeq_16SOTU", package = "DspikeIn")

relative_Des <- physeq_16SOTU %>%
  phyloseq::subset_taxa(Genus != "Tetragenococcus") %>%
  phyloseq::subset_samples(Clade.Order == "Caudate") %>%
  phyloseq::subset_samples(Host.genus %in% c("Desmognathus", "Plethodon", "Eurycea")) %>%
  phyloseq::subset_samples(Ecoregion.III == "Blue Ridge") %>%
  phyloseq::subset_samples(Host.genus == "Desmognathus")


results_DESeq2_rel <- perform_and_visualize_DA(
  obj = relative_Des,
  method = "DESeq2",
  group_var = "Host.taxon",
  contrast = c("Desmognathus monticola", "Desmognathus imitator"),
  output_csv_path = "DA_DESeq2_rel.csv",
  target_glom = "Genus",
  significance_level = 0.05
)

print(results_DESeq2_rel$plot)
head(results_DESeq2_rel$results)
results_DESeq2_rel$obj_significant


```

| Absolute FDR | Relative FDR |
|:---------------------:|:--------------:|
| ![Absolute FDR](https://github.com/user-attachments/assets/0afd03bb-6ee8-450f-aae5-b9146ae522ac) | ![Relative FDR](https://github.com/user-attachments/assets/37c3d8e2-a807-460f-8b6a-9d2c13f0ffc9) |


---

### Customized filtering

```r
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

## Visualization

```r
# taxa barplot
# taxa barplot
 #abundance_type = "absolute"/"relative"
 #ps_ABS= phyloseq object with absolute counts
 
 
 bp_free <- taxa_barplot( physeq = ps_ABS,
                          target_glom = "Genus",
                          treatment_variable = "Genotype",
                          fill_variable = "Genus",           # Fill bars by Genus
                          abundance_type = "absolute",
                          facet_variable = "Treatment",      # Facet by Treatment 
                          x_scale = "free",                  # can be "fixed"
                          legend_size = 10,                
                          top_n_taxa = 30,                   # the top 30 taxa
                          xlab=NULL,                          
                          legend_columns = 1,
                          palette = color_palette$MG)  # This is DspikeIn custom color palette (MG)
 print(bp_free$barplot)
 
 #ps_rel= phyloseq object with relative counts
 # x_scale = "fixed" or "free" 
 
 bp_fix <- taxa_barplot(physeq = ps_rel,
                        target_glom = "Genus",
                        treatment_variable = "Genotype",
                        abundance_type = "relative",
                        facet_variable = "Treatment",
                        fill_variable = "Genus",
                        x_scale = "free", 
                        legend_size = 10,
                        top_n_taxa = 30,
                        xlab=NULL,
                        legend_columns = 1,    # legend col
                        palette = color_palette$MG)
 print(bp_fix$barplot)
 
 
```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Absolute taxa barplot](https://github.com/user-attachments/assets/6d32562f-783e-4462-bd86-5ea1466b5797) | ![Relative taxa barplot](https://github.com/user-attachments/assets/c8a237e9-f707-4dd9-aca6-60cf3c5a6aef) |



```r

# Check abundance distribution via Ridge Plots before and after converting to absolute abundance
ridgeP_before <- ridge_plot_it(spiked_16S_OTU, taxrank = "Family", top_n = 10)
ridgeP_after <- ridge_plot_it(physeq_absolute_16S_OTU, taxrank = "Family", top_n = 10)+ my_custom_theme() #custome theme embedded in DspikeIn

```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/7b50556d-77d7-4aa0-b2b0-f461c67c65a7) | ![Rel  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/5b5a85ea-1f10-4082-a108-7c9019bd84d8) |

```r

# core_microbiome
#ps_ABS=absolute count
#ps_Rel=relative count

custom_detections <- list(
  prevalences = seq(0.05, 1, 0.01),  # Custom prevalences
  thresholds = 10^seq(log10(0.05), log10(1), length = 8),  # detection thresholds
  min_prevalence = 0.4,  # Custom minimum prevalence
  taxa_order = "ascending")  # Order taxa by ascending/desc abundance


# Create and save the plot
plot_result <- plot_core_microbiome_custom(physeq = ps_ABS, 
                                           detections = custom_detections, 
                                           taxrank = "Family", 
                                           output_core_rds = "core_microbiome.rds", 
                                           output_core_csv = "core_microbiome.csv")
print(plot_result)

plot_result <- plot_core_microbiome_custom(physeq = ps_Rel, 
                                           detections = custom_detections, 
                                           taxrank = "Family", 
                                           output_core_rds = "core_microbiome.rds", 
                                           output_core_csv = "core_microbiome.csv")

# core.microbiome is automatically saved in your working directory so you can go ahead and plot it in any way you prefer
core.microbiome <- readRDS("core.microbiome.rds")

```



| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs core](https://github.com/mghotbi/DspikeIn/assets/29090547/9f32f799-4421-4842-a4cf-a1eee1a768e1) | ![Rel core](https://github.com/mghotbi/DspikeIn/assets/29090547/5fb4055a-33f5-48ec-bb0c-a9a9ab446429) |


### Visualization Core microbiome distribution
```r

# shift to long-format data frame and plot the abundance of taxa across the factor of your interest
# Generate alluvial plot
# Convert a phyloseq object to a long-format data frame
# pps_Abs <- phyloseq::psmelt(ps_ABS)

# Example of total reads calculation for relative abundance
# total_reads <- sum(pps_Abs$Abundance)

# Generate an alluvial plot
# Alluvial plot accepts facet

alluvial_plot_abs <- alluvial_plot( data = pps_Abs,
axes = c("Clade.Order","Animal.ecomode","Ecoregion.III", "Diet", "Reproduction","Metamorphosis"),
abundance_threshold = 10000,
fill_variable = "Genus",
silent = TRUE,
abundance_type = "absolute",
top_taxa = 23,
text_size = 4,
legend_ncol = 1,
custom_colors = DspikeIn::color_palette$light_MG)  # Use the extended palette from your package


# Generate alluvial plot for relative abundance
 pps_Rel <- phyloseq::psmelt(ps_Rel)
total_reads <- sum(pps_Rel$Abundance)

alluvial_plot_rel <- alluvial_plot(data = pps_rel,
axes = c("Clade.Order","Animal.ecomode","Ecoregion.III", "Diet", "Reproduction","Metamorphosis"),
abundance_threshold = 0.1,
fill_variable = "Genus",
silent = TRUE,
abundance_type = "relative",
total_reads = total_reads,
top_taxa = 23,
 text_size = 4,
 legend_ncol = 1,
custom_colors = DspikeIn::color_palette$light_MG)


```



| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Absolute](https://github.com/user-attachments/assets/3e7b0ade-c119-4c64-9d96-5638fc7c2296) | ![Relative](https://github.com/user-attachments/assets/a6e75fec-8590-4b58-87f7-aa12a348002a) |

```r

# selecting the most important ASVs/OTUs through RandomForest classification
# Salamander_absolute= subset of our phyloseq object
rf_physeq <- RandomForest_selected(ps_physeq_absolute_16S_OTU, response_var = "Host_Species", na_vars = c("Habitat","Diet", "Ecoregion_III", "Host_genus", "Animal_type"))
RP=ridge_plot_it(rf_physeq)
RP+facet_wrap(~Diet)


```


## Detect common ASVs-OTUs
### DspikeIn includes a feature for detecting common taxa that are dominant or essential across multiple methods.

```r

#detect common ASVs/OTUs
# The input is the list of phyloseq objects
results <- detect_common_asvs_taxa(list(rf_physeq, FTspiked_16S , core.microbiome), 
                                    output_common_asvs_rds = "common_asvs.rds", 
                                    output_common_taxa_rds = "common_taxa.rds")


```

## If you use this package and find it useful:

- ⭐ **Give us stars on GitHub!**
- 🖋️ **Cite us using the link below:**

[Citing DspikeIn](https://www.biorxiv.org/content/10.1101/2024.12.27.630554v1)

![Thank You](https://img.shields.io/badge/Thank-You-brightgreen?style=for-the-badge)

---


