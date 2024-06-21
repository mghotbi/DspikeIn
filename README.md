# DspikeIn
The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package

![How it works](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/Copy%20of%20Untitled%20(2).png)


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
- These DNA isolations served as standards to determine the appropriate spike-in volume of cells to represent 0.1-10% of a sample, as detailed in [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8).

---

## Bioinformatics and Downstream Analyses

### Gene Marker Analysis
- 16S rRNA and ITS rDNA gene markers were analyzed.

### Normalization
- Normalization was performed on the 16S community-weighted mean ribosomal operon copy numbers to correct for potential biases in the representation of *Tetragenococcus halophilus* species in our ASV approach.


````r


---

## Using QIIME2 Plugin for GCN Normalization

To normalize your data by gene copy number (GCN) using the QIIME2 plugin, follow the steps below. For more information, visit the [q2-gcn-norm GitHub repository](https://github.com/Jiung-Wen/q2-gcn-norm).

### Command

Run the following command to perform GCN normalization:

```bash
qiime gcn-norm copy-num-normalize \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --o-gcn-norm-table table-normalized.qza
```

---

This version uses markdown formatting to enhance readability. The link to the GitHub repository is embedded in the text, and the command is formatted as a code block to make it stand out.

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

*The full dataset is available upon request. A subset of the dataset is attached for use in this workshop.*

To test the DspikeIn package, you can download the dataset using the following link: [Download Dataset](https://drive.google.com/file/d/1Ohac-RnrXWSuBAMfqxVzrxCE3Sq98LJK/view?usp=sharing).


```r

# Installation
#Instructions for how to install the DspikeIn package.

install.packages("devtools")
devtools::install_github("mghotbi/DspikeIn")
library(DspikeIn)

# Optional Package Installation
# For users convenience, we provide a helper function to install and load several microbial-ecology-relevant packages, some of which are required for running the `DspikeIn` package.
# You can use the `install_load_packages_helper()` function to easily install and load these packages.
install_load_packages_helper()

```



```r
# Make a new directory and set it as your working directory
create_directory("DspikeIn_16S_OTU", set_working_dir = TRUE)
getwd()

# Please note that these functions have been primarily built using several wonderful packages, including phyloseq (https://github.com/joey711/phyloseq) and tidyverse/dplyr (https://github.com/tidyverse/dplyr).
# Therefore, please start by creating a phyloseq object and follow the instructions.
# To create your phyloseq object, please refer to the phyloseq tutorial (https://joey711.github.io/phyloseq).
# The phyloseq object needs to include OTU/ASV, Taxa, phylogenetic tree, DNA reference, 
# and metadata containing spiked species volume, starting from 0 (no spike species added) to 4 (4 μl of spike cell added).

# Note: DspikeIn requires 'spiked.volume'; any other format is not readable.
print_sentence("¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯")


# We are going to work with a subset of the dataset for both ASVs and OTUs approaches to accelerate this workshop.

Salamander_relative_16S_ASV <-readRDS("Salamander_relative_16S_ASV.rds")
Salamander_relative_ITS_ASV <-readRDS("Salamander_relative_ITS_ASV.rds")

Salamander_relative_16S_OTU <-readRDS("Salamander_relative_16S_OTU.rds")
Salamander_relative_ITS_OTU <-readRDS("Salamander_relative_ITS_OTU.rds")


physeq_16S_ASV <- tidy_phyloseq(Salamander_relative_16S_ASV)

# Ensure your metadata contains spiked volumes:
physeq_16S_ASV@sam_data$spiked.volume


```


# Prepare the required information 


```r
# Required Information 
# Please note that the Spike cell numbers, species name, and selected hashcodes are customizable and can be tailored to the specific needs of individual studies.
# Moreover, to proceed with the DspikeIn package, you only need to select one method to specify your spiked species: either by hashcodes or species name.

# 16S rRNA
spiked_cells <-1847
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
merged_spiked_species<-"Tetragenococcus_halophilus"
Tetra <- subset_taxa(physeq_16SASV,Species=="Tetragenococcus_halophilus" | Species=="Tetragenococcus_sp")
hashcodes <- row.names(phyloseq::tax_table(Tetra))

# ITS rDNA
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
# using the Neighbor-Joining method  based on a Jukes-Cantor distance matrix and plot the tree with bootstrap values.
# we compare the Sanger read of Tetragenococcus halophilus with the FASTA sequence of Tetragenococcus halophilus from our phyloseq object.

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
| ![Neighbor Joining Tree](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/JU.png?raw=true) | ![Tetra Plot](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/tetra.png?raw=true) | ![cophenetic](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/Rplot02.png?raw=true) |




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



![WhyOTUs](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/ASVOTU.png)

---

Retrieved spiked species are correlated with abundance and can therefore vary per system of study. The retrieved spiked species across taxa abundance was tested. We divided this range into *0-10%, 10-15%, 15-25%, 25-35%, and 35-100%* to test which range can lead to significant variation in the percentage of retrieved spiked species and scaling factors across abundance. The results revealed that we can expand the acceptable range of retrieved spiked species percentage to 35% in our model system, which contrasts with [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8), who noted that the acceptable range of retrieved spiked species can be between 0.1% and 10%.


---




| ASVs Or OTUs | Acceptable range|
|:----------:|:---------:|
| ![ASVs vs OTUs](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/pngof%20desired.png) | ![Acceptable range](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/range.png) |



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


### Data Transformation

*Check if the transformations is required.*

```r

# Adjust abundance by one-third
readAdj16S <- adjust_abundance_one_third(spiked_16S_OTU, factor = 3)
summ_count_phyloseq(readAdj16S)

# Random subsampling with reduction factor
red16S <- random_subsample_WithReductionFactor(spiked_16S_OTU, reduction_factor = 3)
summ_count_phyloseq(red16S)

# Proportion adjustment
normalized_16S <- proportion_adj(spiked_16S_OTU, output_file = "proportion_adjusted_physeq.rds")
summ_count_phyloseq(normalized_16S)

# DESeq2 variance stabilizing transformation (VST)
transformed_16S <- run_vst_analysis(spiked_16S_OTU)
summ_count_phyloseq(transformed_16S)

# Relativize and filter taxa based on selected thresholds
FTspiked_16S <- relativized_filtered_taxa(
  spiked_16S_OTU,
  threshold_percentage = 0.0001,
  threshold_mean_abundance = 0.0001,
  threshold_count = 5,
  threshold_relative_abundance = 0.0001)
summ_count_phyloseq(FTspiked_16S)

# Adjust prevalence based on the minimum reads
spiked_16S_min <- adjusted_prevalence(spiked_16S_OTU, method = "min")


```


## Preprocessing for Scaling Factor Calculation  

If you are using OTUs and have only one OTU rooted from the spiked species, you can skip this preprocessing step. Follow the steps below to estimate the success of spike-in, particularly check if you have any samples with under or over-spikes.If the spiked species appear in several ASVs, check their phylogenetic distances and compare them to the reference sequences of your positive control. If the spiked species of interest has gene copy number variations and you prefer not to sum their abundances, use the `max` option instead of `sum` to combine these ASVs under a single taxon, simplifying data processing.

```


# Modify the threshold of acceptable spiked species % as needed. For detailed guidance on acceptable thresholds (passed_range), please refer to the instructions in our upcoming paper.

# Merg the spiked species
# merge_method = "max": Selects the maximum abundance among ASVs of the spiked species, ensuring the most abundant ASV is retained.
# merge_method = "sum": Sums the abundances of ASVs of the spiked species, providing a cumulative total.


species_name <- "Tetragenococcus_halophilus"
Spiked_16S_sum_scaled <- Pre_processing_species(spiked_16S_OTU, species_name, merge_method = "sum", output_file = "merged_physeq_sum.rds")
Spiked_16S_max_scaled <- Pre_processing_species(spiked_16S_OTU, species_name, merge_method = "max", output_file = "merged_physeq_max.rds")

Spiked_16S_sum_scaled <- Pre_processing_hashcodes(spiked_16S_OTU, hashcodes, merge_method = "sum", output_prefix = "merged_physeq_sum")
Spiked_16S_max_scaled <- Pre_processing_hashcodes(spiked_16S_OTU, hashcodes, merge_method = "max", output_prefix = "merged_physeq_max")
summ_count_phyloseq(Spiked_16S_sum_scaled)


Spiked_16S_OTU_scaled <- tidy_phyloseq(Spiked_16S_sum_scaled)
# Now calculate the spiked species retrieval percentage. Customize the passed_range and merged_spiked_species/merged_spiked_hashcodes based on your preferences.
# passed_range = "c(0.1, 10) ": threshold of acceptable spiked species %
# passed_range = "c(0.1, 35) ": threshold of acceptable spiked species %
## select either merged_spiked_species or merged_spiked_hashcodes
# merged_spiked_species = merged_spiked_species 
# merged_spiked_hashcodes = merged_spiked_hashcodes

merged_spiked_species <- c("Tetragenococcus_halophilus")
result <- calculate_spike_percentage(Spiked_16S_OTU_scaled, merged_spiked_species, passed_range = c(0.1, 11))
calculate_summary_stats_table(result)

# Define your merged_spiked_hashcodes
merged_Tetra <- subset_taxa(Spiked_16S_OTU_scaled, Species == "Tetragenococcus_halophilus")
merged_spiked_hashcodes <- row.names(tax_table(merged_Tetra))
result <- calculate_spike_percentage(Spiked_16S_OTU_scaled,  merged_spiked_hashcodes, passed_range = c(0.1, 35))
calculate_summary_stats_table(result)



# If you decide to remove the failed reads and go forward with passed reads, here is what you need to do
# you can also go forward with the original file and remove the failed reads after converting relative to absolute abundance
# Filter to get only the samples that passed
passed_samples <- result$Sample[result$Result == "passed"]
# Subset the original phyloseq object to keep only the samples that passed
passed_physeq <- prune_samples(passed_samples, Spiked_16S_ASV_scaled)

```

## Data Normalization and Transformation
*Experiment Repetition*

Getting help from [Yerk et al., 2024](https://doi.org/10.1186/s40168-023-01747-z), we checked if we needed to normalize our data before/after calculating our spiked species to account for spiked volume variations and library preparation. We evaluated the need for compositionally aware data transformations, including centered log-ratio (CLR) transformation, DESeq2 variance stabilizing transformation (`run_vst_analysis`), subsampling with a reduced factor for count data (`random_subsample_WithReductionFactor`), proportion adjustment (`proportion.adj`), and prevalence adjustment (`adjusted_prevalence`). Additionally, we considered compositionally naïve data transformations, such as raw data and relative abundance-based transformations (`relativized_filtered_taxa`), before calculating spike-in scaling factors. The only significant variation in the percentage of retrieved spiked species was relevant to VST, so we continued with raw data.


You can repeat the experiment by transforming the data, calculating spike percentage using `calculate_spike_percentage()`, then checking for the homogeneity of variances using `Bartlett_test()` and ensuring the data is normally distributed using `Shapiro_Wilk_test()`. Finally, plot the results using `transform_plot()`.

```r

methods<-readRDS("methods.rds")
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

# y_vars are numerical variables of your interest to be analysed
y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")
# x_var is a categorical variable
x_var <- "Methods"
# the color_pallet is MG here

# scale data
scaled <- methods %>% mutate_at(c("Total.reads", "Spike.reads", "Spike.percentage" ), ~(scale(.) %>% as.vector))

# Perform Kruskal-Wallis test
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", color_palette = MG, stat_test = "anova")
# Perform one-way ANOVA
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", color_palette = MG, stat_test = "kruskal.test")


```


| Spiked sp Percentage ANOVA | Spiked sp Reads ANOVA | Total Reads ANOVA |
|:--------------------------:|:---------------------:|:-----------------:|
| ![Spike Percentage ANOVA](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/plot_Spike.percentage_ANOVA.jpg?raw=true) | ![Spike Reads ANOVA](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/plot_Spike.reads_ANOVA.png?raw=true) | ![Total Reads ANOVA](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/plot_Total.reads_ANOVA.jpg?raw=true) |

---



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

# summary statistics 
post_eval_summary <- calculate_summary_stats_table(physeq_16S_adj_scaled_AbsoluteCount)
print(post_eval_summary)

# Create a new phyloseq obj with absolute counts
physeq_16S_adj_scaled_absolute_abundance <- phyloseq(
  otu_table = round(otu_table(Spiked_16S_OTU_scaled) * scaling_factors), 
  taxa_table = tax_table(Spiked_16S_OTU_scaled),
  phy_tree = phy_tree(Spiked_16S_OTU_scaled),
  sample_data = sample_data(Spiked_16S_OTU_scaled))

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
![Success](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/image.png)


*Save your file for later. Stay tuned for the rest: transformation and several visualization methods and displaying the importance of converting relative to absolute abundance in the context of microbial ecology.*


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
| ![Absolute Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/Sal.abs.pnnn.png) | ![Relative Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/sal.rel.pmnb.png) |




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
| ![Absolute Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/ridge%20abs1.png) | ![Relative Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/ridge%20rel2.png) |


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
| ![Absolute Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/Abs.coree.png) | ![Relative Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/core.rel.png) |


```r

# shift to dataframe and plot the abundance of taxa across the factors
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
| ![Absolute Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/abs.eq0.png) | ![Relative Abundance](https://github.com/mghotbi/DspikeIn/blob/MGhotbi/eq%20rel3.png) |


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

plot_asvs_abundance(common_asvs_phyloseq, response_var = "host.species", x_var = "ecoregion.III", rank_var = "Phylum")

```
