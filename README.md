# DspikeIn
The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package

![DspikeIn](https://github.com/mghotbi/DspikeIn/assets/29090547/f827aed9-2f99-42a1-adff-4b85124f8e94)

---

## Introduction

In our study, *Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples based on our previous studies [WalkerLab](https://walkerlabmtsu.weebly.com/personnel.html).



## Methodology

### Growth of Stock Cell Suspensions
- _**Tetragenococcus halophilus**_: Cultivated in tryptic soy broth.
- _**Dekkera bruxellensis**_: Cultivated in potato dextrose broth.
- Both microbial cultures were serially diluted, and optical density (OD) measurements were obtained using a ClarioStar plate reader. To obtain an aliquot and follow our procedure, please contact [Prof. Donnald Walker](mailto:Donald.Walker@mtsu.edu).

### DNA Extraction

- DNA was extracted using the Qiagen DNeasy Powersoil Pro Kit.
- These DNA isolations served as standards to determine the appropriate spike-in volume of cells to represent 0.1-10% of a sample, as detailed in [Rao et al., 2021](https://www.nature.com/articles/s41586-021-03241-8).

---

## Bioinformatics and Downstream Analyses

### Gene Marker Analysis
- 16S rRNA and ITS rDNA gene markers were analyzed.

### Normalization
- Normalization was performed on the 16S community-weighted mean ribosomal operon copy numbers to correct for potential biases in the representation of *Tetragenococcus halophilus* species in our ASV approach.


```markdown

## Using QIIME2 Plugin for GCN Normalization

# To normalize data by gene copy number (GCN) using the QIIME2 plugin, follow the steps below.
# For more information, visit the q2-gcn-norm GitHub repository (https://github.com/Jiung-Wen/q2-gcn-norm).

### Command

Run the following command to perform GCN normalization:

qiime gcn-norm copy-num-normalize \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --o-gcn-norm-table table-normalized.qza

```

---


### DspikeIn Package
The **DspikeIn** package was developed to facilitate:
- Verifying the phylogenetic distances of ASVs/OTUs rooted from spiked species.
- Preprocessing data.
- Calculating the spike-in scaling factor.
- Converting relative abundance to absolute abundance.
- Data transformation and visualization.

We provide a step-by-step walkthrough of each procedure within the **DspikeIn** package.

For more detailed methodology and results, please refer to our soon-to-be-published paper.

---



Similar to the debate on the application of OTUs vs. ASVs, which has many contrasting views on their benefits and drawbacks ([Callahan et al., 2017](https://doi.org/10.1038/ismej.2017.119); [Schloss et al., 2021](https://doi.org/10.1128/msphere.00191-21)), there are positive ([Schirrmeister et al., 2012](https://doi.org/10.1186/1471-2180-12-177); [Stoddard et al., 2015](https://doi.org/10.1093/nar/gku1201)) and opposing ([Louca et al., 2018](https://doi.org/10.1038/s43705-023-00266-0); [Gao and Wu et al., 2023](https://doi.org/10.1038/s43705-023-00266-0)) opinions on gene copy number correction for the 16S rRNA marker. Meanwhile, several novelties and modifications have been added to copy number correction to improve accuracy ([Perisin et al., 2016](https://doi.org/10.1038/ismej.2015.161); [Gao and Wu, 2023](https://doi.org/10.1038/s43705-023-00266-0)).



Using our spike-in positive controls and assessing the percentage of retrieved spiked species, we were able to compare the results at each step before selecting a pathway. The copy number correction was performed to correct potential biases in the representation of *Tetragenococcus halophilus* species in our ASV approach. We utilized `q2-gcn-norm` based on the rrnDB database (version 5.7) to normalize for 16S rRNA gene marker copy numbers ([qiime2 plugin; gcn-norm](https://github.com/Jiung-Wen/q2-gcn-norm)). Due to the variability in rDNA gene copy numbers, straightforward translation of rDNA read counts into the abundance of individual organisms is precluded ([Lavrinienko et al., 2021](https://doi.org/10.1016/j.tim.2020.05.019)).



In fact, for ITS we did not need to use copy number correction. However, we recommend conducting a literature review before deciding whether to sum or select the maximum abundance of OTUs/ASVs rooted from spiked-in species and calculating spike-in factors. We believe systematic evaluation before selecting or refuting each method can help prevent miscalculations ([Lofgren et al., 2018](https://doi.org/10.1111/mec.14995)), which can aid in establishing a system-dependent method for copy number correction in ITS markers.

---
## Dataset

*The full dataset will be available upon request. A subset of the dataset is attached for use in this workshop.*

To test the DspikeIn package, you can download the dataset using the following link: [Download Dataset](https://drive.google.com/file/d/1Ohac-RnrXWSuBAMfqxVzrxCE3Sq98LJK/view?usp=sharing).



*If you encounter issues installing the package due to missing dependencies, follow these steps to install all required packages first:*

## Step 1: Install Required Packages

To install the required packages, use the following script:

```r

# Please install devtools first
install.packages("devtools")

# Source the helper function from our GitHub repository
devtools::source_url("https://raw.githubusercontent.com/mghotbi/DspikeIn/MGhotbi/install_required_packages.R")

# Run this function to install all required packages
install_required_packages()


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


## Optional Package Installation
# For users convenience, we provide a helper function to install and load several microbial-ecology-relevant packages, some of which are required for running the `DspikeIn` package.
# You can use the `install_load_packages_helper()` function to easily install and load these packages.
install_load_packages_helper()

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

# Note: DspikeIn requires 'spiked.volume'; any other format is not readable.
print_sentence("¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯")


# We are going to work with a subset of the dataset for both ASVs and OTUs
# approaches to accelerate this workshop.

Salamander_relative_16S_ASV <-readRDS("Salamander_relative_16S_ASV.rds")
Salamander_relative_ITS_ASV <-readRDS("Salamander_relative_ITS_ASV.rds")

physeq_16S_ASV <- tidy_phyloseq(Salamander_relative_16S_ASV)

# Ensure your metadata contains spiked volumes:
physeq_16S_ASV@sam_data$spiked.volume


```


# Prepare the required information 


```r
# Required Information 
# Please note that the Spike cell numbers, species name, and selected hashcodes are customizable and can be tailored to the specific needs of individual studies.
# Moreover, to proceed with the DspikeIn package, you only need to select one method to specify your spiked species: either by hashcodes or species name.

library(phyloseq)
# 16S rRNA
presence of 'spiked.volume' column in metadata
spiked_cells <-1847
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
merged_spiked_species<-"Tetragenococcus_halophilus"
Tetra <- subset_taxa(physeq_16SASV,Species=="Tetragenococcus_halophilus" | Species=="Tetragenococcus_sp")
hashcodes <- row.names(phyloseq::tax_table(Tetra))

# ITS rDNA
presence of 'spiked.volume' column in metadata
spiked_cells <- 733
species_name <- spiked_species<-merged_spiked_species<-"Dekkera_bruxellensis"
Dekkera <- subset_taxa(physeq_ITSASV, Species=="Dekkera_bruxellensis")
hashcodes <- row.names(phyloseq::tax_table(Dekkera))

```


# Plot phylogenetic tree with Bootstrap Values
This step will be helpful for handling ASVs with/without Gene Copy Number Correction
This section demonstrates how to use various functions from the package to plot and analyze phylogenetic trees.



```r
# In case there are still several ASVs rooting from the spiked species, you may want to check the phylogenetic distances.
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



![WhyOTUvsASV](https://github.com/mghotbi/DspikeIn/assets/29090547/10d82aea-9aa5-476e-a421-e0e6ddb89841)


---

Retrieved spiked species are correlated with abundance and can therefore vary per system of study. The retrieved spiked species across taxa abundance was tested. We divided this range into *0-10%, 10-15%, 15-25%, 25-35%, and 35-100%* to test which range can lead to significant variation in the percentage of retrieved spiked species and scaling factors across abundance. The results revealed that we can expand the acceptable range of retrieved spiked species percentage to 35% in our model system, which contrasts with [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8), who noted that the acceptable range of retrieved spiked species can be between 0.1% and 10%.


---




| ASVs Or OTUs | Acceptable range|
|:----------:|:---------:|
| ![Desired range of spiked sp](https://github.com/mghotbi/DspikeIn/assets/29090547/2f949616-6493-4445-8f1e-7ac9c9dd844f) | ![Acceptable range](https://github.com/mghotbi/DspikeIn/assets/29090547/8674f3de-ba24-4857-9cd7-d1b6dc15c669) |



# By now, we have an idea of what works well for our approach, whether it's OTUs or ASVs.

We will proceed with OTUs and a 35% acceptable range of spiked species retrieval. To learn more about why we chose OTUs over ASVs, and what percentage of retrieved spiked species can be considered passed or failed and why, please read our soon-to-be-published paper. For now, check the schematic of our experiment below.


_The VSEARCH with de novo robust clustering algorithms at a 97% similarity threshold was used to reduce potential mistakes_ [Westcott and Schloss 2015](https://doi.org/10.7717/peerj.1487).



## Subsetting and Preprocessing Spiked Data


*Subset the part of the data which is spiked. Keep solely spiked samples using the `spiked.volume` column.*


```r

# Subset spiked samples (264 samples are spiked)
spiked_16S_OTU <- subset_samples(physeq_16S_OTU, spiked.volume %in% c("2", "1"))
spiked_16S_OTU <- tidy_phyloseq(spiked_16S_OTU)

```

### Examine Your Count Table/Biom File Before Going Further


```r

# Summarize the initial statistics for ASVs/OTUs
initial_stat_ASV <- summ_phyloseq_ASV_OTUID(physeq_16S_OTU)

# Summarize the initial statistics sample-wise
initial_stat_sampleWise <- summ_phyloseq_sampleID(physeq_16S_OTU)

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

If you are using OTUs and have only one OTU rooted from the spiked species, you can skip this preprocessing step. Follow the steps below to estimate the success of spike-in, particularly check if you have any samples with under or over-spikes.If the spiked species appear in several ASVs, check their phylogenetic distances and compare them to the reference sequences of your positive control. If the spiked species of interest has gene copy number variations and you prefer not to sum their abundances, use the `max` option instead of `sum` to combine these ASVs under a single taxon, simplifying data processing.


```r

# Modify the threshold of acceptable spiked species % as needed. 
# For detailed guidance on acceptable thresholds (passed_range), 
# please refer to the instructions in our upcoming paper.

# Merge the spiked species
# merge_method = "max": Selects the maximum abundance among ASVs of the spiked species, 
# ensuring the most abundant ASV is retained.
# merge_method = "sum": Sums the abundances of ASVs of the spiked species, 
# providing a cumulative total.

species_name <- "Tetragenococcus_halophilus"

# Merge using "sum" method
Spiked_16S_sum_scaled <- Pre_processing_species(
  spiked_16S_OTU, 
  species_name, 
  merge_method = "sum", 
  output_file = "merged_physeq_sum.rds")

# Merge using "max" method
Spiked_16S_max_scaled <- Pre_processing_species(
  spiked_16S_OTU, 
  species_name, 
  merge_method = "max", 
  output_file = "merged_physeq_max.rds")

# Merge hashcodes using "sum" method
Spiked_16S_sum_scaled <- Pre_processing_hashcodes(
  spiked_16S_OTU, 
  hashcodes, 
  merge_method = "sum", 
  output_prefix = "merged_physeq_sum")

# Merge hashcodes using "max" method
Spiked_16S_max_scaled <- Pre_processing_hashcodes(
  spiked_16S_OTU, 
  hashcodes, 
  merge_method = "max", 
  output_prefix = "merged_physeq_max")

# Summarize count
summ_count_phyloseq(Spiked_16S_sum_scaled)

# Tidy phyloseq object
Spiked_16S_OTU_scaled <- tidy_phyloseq(Spiked_16S_sum_scaled)

# Now calculate the spiked species retrieval percentage.
# Customize the passed_range and merged_spiked_species/merged_spiked_hashcodes based on your preferences.
# passed_range = "c(0.1, 10)": threshold of acceptable spiked species %
# passed_range = "c(0.1, 35)": threshold of acceptable spiked species %
# Select either merged_spiked_species or merged_spiked_hashcodes

merged_spiked_species <- c("Tetragenococcus_halophilus")
result <- calculate_spike_percentage(
  Spiked_16S_OTU_scaled, 
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
  passed_range = c(0.1, 35))
calculate_summary_stats_table(result)

# If you decide to remove the failed reads and go forward with passed reads, here is what you need to do
# You can also go forward with the original file and remove the failed reads 
# after converting relative to absolute abundance

# Filter to get only the samples that passed
passed_samples <- result$Sample[result$Result == "passed"]

# Subset the original phyloseq object to keep only the samples that passed
passed_physeq <- prune_samples(
  passed_samples, 
  Spiked_16S_ASV_scaled)

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


# Convert Relative Counts to Absolute Counts and Create a New Phyloseq Object


```r

# Convert relative counts data to absolute counts
physeq_16S_adj_scaled_AbsoluteCount <- convert_to_absolute_counts(Spiked_16S_OTU_scaled, scaling_factors)
absolute <- convert_to_absolute_counts(Spiked_16S_OTU_scaled, scaling_factors)
absolute_counts <- physeq_16S_adj_scaled_AbsoluteCount$absolute_counts
physeq_absolute_abundance_16S_OTU <- physeq_16S_adj_scaled_AbsoluteCount$physeq_obj


# summary statistics 
post_eval_summary <- calculate_summary_stats_table(absolute_counts)
print(post_eval_summary)


```


# Let's check the conclusion and get the report table of spiked species success or failure.



```r

# Define the parameters once.
merged_spiked_species <- c("Tetragenococcus_halophilus")
max_passed_range <- 35

# Subset the phyloseq object to exclude blanks
physeq_16S_adj_scaled_perc <- subset_samples(Spiked_16S_OTU_scaled, sample.or.blank != "blank")

# Generate the spike success report and summary statistics
summary_stats <- conclusion(physeq_16S_adj_scaled_perc, merged_spiked_species, max_passed_range)
print(summary_stats)


```

Here is an example of a success or failure report:
![success report](https://github.com/mghotbi/DspikeIn/assets/29090547/017cfa65-8b75-4625-8d49-6e4a67146193)



```r

#Save your file for later. Please stay tuned for the rest: Comparisons and several visualization methods to show how important it is to convert relative to absolute abundance in the context of microbial ecology.

taxa_names(physeq_absolute_abundance_16S_OTU) <- paste0("ASV", seq(ntaxa(physeq_absolute_abundance_16S_OTU)))
physeq_absolute_abundance_16S_OTU <- tidy_phyloseq(physeq_absolute_abundance_16S_OTU)
saveRDS(physeq_absolute_abundance_16S_OTU, "physeq_absolute_abundance_16S_OTU.rds")

```
# Normalization and bias correction 


```r
# Bolstad, B.M., Irizarry, R.A., Åstrand, M. and Speed, T.P., 2003. A comparison of normalization methods for high-density oligonucleotide array data based on variance and bias. Bioinformatics, 19(2), pp.185-193.
# Gagnon-Bartsch, J.A. and Speed, T.P., 2012. Using control genes to correct for unwanted variation in microarray data. Biostatistics, 13(3), pp.539-552.
# Risso, D., Ngai, J., Speed, T.P. and Dudoit, S., 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature biotechnology, 32(9), pp.896-902.
# Gagnon-Bartsch, J.A., Jacob, L. and Speed, T.P., 2013. Removing unwanted variation from high dimensional data with negative controls. Berkeley: Tech Reports from Dep Stat Univ California, pp.1-112.

# Load required libraries
library(phyloseq)
library(DESeq2)
library(edgeR)
library(preprocessCore)
library(sva)
library(EDASeq)
library(Biobase)
library(BiocGenerics)
library(vegan)


#ps is a phyloseq object without spiked species counts
#One can calculate the scaling factor using any normalization method in the absence of spiked species counts, and then determine the spiked scaling factor. Crossing both #scaling factors with relative abundance helps quantify absolute abundance while correcting for bias

ps <- remove_zero_negative_count_samples(physeq_absolute_abundance_16S_OTU)
ps <- convert_categorical_to_factors(physeq_absolute_abundance_16S_OTU)

# Normalization Methods:
# group_var <- "Animal.ecomode"  
# result_DESeq <- normalization_set(ps, method = "DESeq", groups = group_var)
# result_TMM <- normalization_set(ps, method = "TMM", groups = group_var)
# result_CLR <- normalization_set(ps, method = "clr")
# result_SVA <- normalization_set(ps, method = "SVA", groups = group_var)
# result_RUVg <- normalization_set(ps, method = "RUVg", groups = group_var)
# result_RUVr <- normalization_set(ps, method = "RUVr", groups = group_var)
# result_RUVs <- normalization_set(ps, method = "RUVs", groups = group_var)
# result_UQ <- normalization_set(ps, method = "UQ", groups = group_var)
# result_med <- normalization_set(ps, method = "med", groups = group_var)
# result_rle <- normalization_set(ps, method = "rle")
# result_css <- normalization_set(ps, method = "CSS")
# result_tss <- normalization_set(ps, method = "tss")
# result_rar <- normalization_set(ps, method = "rar")

# Customized filtering and transformations
# Proportion adjustment
physeq<-physeq_absolute_abundance_16S_OTU
normalized_physeq <- proportion_adj(physeq, output_file = "proportion_adjusted_physeq.rds")
summ_count_phyloseq(normalized_16S)


# Relativize and filter taxa based on selected thresholds
FT_physeq <- relativized_filtered_taxa(
  physeq,
  threshold_percentage = 0.0001,
  threshold_mean_abundance = 0.0001,
  threshold_count = 5,
  threshold_relative_abundance = 0.0001)
summ_count_phyloseq(FT_physeq)

# Adjust prevalence based on the minimum reads
physeq_min <- adjusted_prevalence(physeq, method = "min")


```

*Experiment Repetition*

Getting help from [Yerk et al., 2024](https://doi.org/10.1186/s40168-023-01747-z), We evaluated the need for compositionally aware data transformations, including centered log-ratio (CLR), and additive log-ratio (alr) transformation, DESeq2 variance stabilizing transformation (`run_vst_analysis`), subsampling with a reduced factor for count data (`random_subsample_WithReductionFactor`), proportion adjustment (`proportion.adj`), and prevalence adjustment (`adjusted_prevalence`). Additionally, we considered compositionally naïve data transformations, such as raw data and relative abundance-based transformations (`relativized_filtered_taxa`), and compared the results. The only noticeable variation in the percentage of retrieved spiked species was related to VST. However, this variation was not significant for spiked sp reterival%.


You can repeat the experiment by transforming the data, calculating spike percentage using `calculate_spike_percentage()`, then checking for the homogeneity of variances using `Bartlett_test()` and ensuring the data is normally distributed using `Shapiro_Wilk_test()`. Finally, plot the results using `transform_plot()`.


```r

methods <- readRDS("methods.rds")
methods$Total.reads <- as.numeric(gsub(",", "", methods$Total.reads))
methods$Spike.reads <- as.numeric(gsub(",", "", methods$Spike.reads))

# Ensure grouping variable is a factor
methods$Methods <- as.factor(methods$Methods)
methods$Result <- as.factor(methods$Result)

# Perform Bartlett test/homogeneity of variances
Bartlett_test(methods, "Result")
Bartlett_test(methods, "Methods")

# Check if data is normally distributed
Shapiro_Wilk_test(methods, "Methods")
Shapiro_Wilk_test(methods, "Result")

# y_vars are numerical variables of your interest to be analyzed
y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")
# x_var is a categorical variable
x_var <- "Methods"
# the color_palette is MG here

# Scale data
scaled <- methods %>% mutate_at(c("Total.reads", "Spike.reads", "Spike.percentage"), ~(scale(.) %>% as.vector))

# Perform Kruskal-Wallis test
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", color_palette = MG, stat_test = "anova")
# Perform one-way ANOVA
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", color_palette = MG, stat_test = "kruskal.test")

```





| Spiked sp Percentage ANOVA | Spiked sp Reads ANOVA | Total Reads ANOVA |
|:--------------------------:|:---------------------:|:-----------------:|
| ![plot_Spike percentage_ANOVA](https://github.com/mghotbi/DspikeIn/assets/29090547/94b4da0b-4dd9-4af9-b897-4207ec2cef46) | ![plot_Spike reads_ANOVA](https://github.com/mghotbi/DspikeIn/assets/29090547/92be2eb3-68c4-4bde-87a9-31821e62c558) | ![plot_Total reads_ANOVA](https://github.com/mghotbi/DspikeIn/assets/29090547/43f7a692-cb15-42c0-b116-f1397619f32d) |

---

## Visualization and Differential abundance 


```r

taxa_names(physeq_16S_adj_scaled_absolute_abundance) <- paste0("ASV", seq(ntaxa(physeq_16S_adj_scaled_absolute_abundance)))
physeq_16S_adj_scaled_absolute_abundance <- tidy_phyloseq(physeq_16S_adj_scaled_absolute_abundance)
saveRDS(physeq_16S_adj_scaled_absolute_abundance, "physeq_16S_adj_scaled_absolute_abundance.rds")


# A subset of the dataset for both relative and absolute abundance of spiked species was filtered

# taxa barplot 
bp_ab <- taxa_barplot(Salamander_absolute_NospikeSp, target_glom = "Genus", treatment_variable = "Host.genus", abundance_type = "absolute", x_angle = 90, fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
print(bp_ab$barplot)

bp_rel <- taxa_barplot(Salamander_relative_NospikeSp, target_glom = "Genus", treatment_variable = "Host.genus", abundance_type = "relative", x_angle = 90, fill_variable = "Genus", facet_variable = "Diet", top_n_taxa = 20)
print(bp_rel$barplot)

```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Rel AbsSal taxa barplot](https://github.com/mghotbi/DspikeIn/assets/29090547/643fca2a-6087-49f1-b3ee-5759d2fcb36f) | ![Rel abun Sal taxa barplot](https://github.com/mghotbi/DspikeIn/assets/29090547/2040830e-1ce1-46c9-8dfe-3c18f77a85bf) |




```r

# simple barplot of taxonomy abundance
# Plot relativized abundance
plot <- plotbar_abundance(Salamander_relative_NospikeSp, level = "Family", group = "Env.broad.scale.x", top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = TRUE, output_prefix = "relativized_abundance_plot")
print(plot)

# Plot non-relativized (absolute) abundance
plot_absolute <- plotbar_abundance(Salamander_absolute_NospikeSp, level = "Family", group = "Env.broad.scale.x", top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE, output_prefix = "non_relativized_abundance_plot")
print(plot_absolute)


# Check abundance distribution via Ridge Plots before and after converting to absolute abundance
ridgeP_before <- ridge_plot_it(Salamander_relative_NospikeSp, taxrank = "Family", top_n = 10)
ridgeP_after <- ridge_plot_it(Salamander_absolute_NospikeSp, taxrank = "Family", top_n = 10)


```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/7b50556d-77d7-4aa0-b2b0-f461c67c65a7) | ![Rel  ridge](https://github.com/mghotbi/DspikeIn/assets/29090547/5b5a85ea-1f10-4082-a108-7c9019bd84d8) |


```r

# core_microbiome
custom_detections <- 10^seq(log10(3e-1), log10(0.5), length = 5)
PCM_rel <- plot_core_microbiome_custom(Salamander_relative_NospikeSp, detections = custom_detections, taxrank = "Family", output_core_rds = "core_microbiome.rds", output_core_csv = "core_microbiome.csv")

PCM_Abs <- plot_core_microbiome_custom(Salamander_absolute_NospikeSp, detections = custom_detections, taxrank = "Family", output_core_rds = "core_microbiome.rds", output_core_csv = "core_microbiome.csv")

# core.microbiome is automatically saved in your working directory so yoou can go ahead and barplot it
core.microbiome <- readRDS("core.microbiome.rds")

```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs core](https://github.com/mghotbi/DspikeIn/assets/29090547/9f32f799-4421-4842-a4cf-a1eee1a768e1) | ![Rel core](https://github.com/mghotbi/DspikeIn/assets/29090547/5fb4055a-33f5-48ec-bb0c-a9a9ab446429) |


```r

# shift to dataframe and plot the abundance of taxa across the factor of your interest
# Load data
meli_Abs_WSal <- readRDS("meli_Abs_WSal.rds")
meli_Rel_WSal <- readRDS("meli_Rel_WSal.rds")

# Generate alluvial plot
is_alluvia_form(as.data.frame(meli_Abs_WSal), axes = 1:6, silent = TRUE)
alluvial_plot <- alluvial_plot(data = meli_Abs_WSal, axes = c("Diet", "Host.genus", "Ecoregion.III"), abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE, abundance_type = "absolute")
print(alluvial_plot)
is_alluvia_form(as.data.frame(meli_Rel_WSal), axes = 1:6, silent = TRUE)
alluvial_plot <- alluvial_plot(data = meli_Rel_WSal, axes = c("Diet", "Host.genus", "Ecoregion.III"), abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE, abundance_type = "relative")
print(alluvial_plot)


```


| Absolute Abundance | Relative Abundance |
|:----------:|:---------:|
| ![Abs Alluv](https://github.com/mghotbi/DspikeIn/assets/29090547/2f187727-db7b-41a2-82be-73162423ce25) | ![Rel Alluv](https://github.com/mghotbi/DspikeIn/assets/29090547/bc6ed255-97d3-4e24-ad22-12890b747e79) |


```r

# selecting the most important ASVs/OTUs through RandomForest classification
# Salamander_absolute= subset of our phyloseq object
rf_physeq <- RandomForest_selected_ASVs(Salamander_absolute, response_var = "Host_Species", na_vars = c("Habitat","Diet", "Ecoregion_III", "Host_genus", "Animal_type"))
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
