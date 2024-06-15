# DspikeIn
The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package

![How it works](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/DspikeIn.png)


---

## Introduction

In our study, *Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples based on our previous studies [WalkerLab](https://walkerlabmtsu.weebly.com/personnel.html).



## Methodology

### Growth of Stock Cell Suspensions
- **Tetragenococcus halophilus**: Cultivated in tryptic soy broth.
- **Dekkera bruxellensis**: Cultivated in potato dextrose broth.
- Both microbial cultures were serially diluted, and optical density (OD) measurements were obtained using a ClarioStar plate reader.

### DNA Extraction

- DNA was extracted using the Qiagen DNeasy Powersoil Pro Kit.
- These DNA isolations served as standards to determine the appropriate spike-in volume of cells to represent 0.1-10% of a sample, as detailed in [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8).

---

## Bioinformatics and Downstream Analyses

### Gene Marker Analysis
- 16S rRNA and ITS rDNA gene markers were analyzed.

### Normalization
- Normalization was performed on the 16S community-weighted mean ribosomal operon copy numbers to correct for potential biases in the representation of *Tetragenococcus halophilus* species in our ASV approach.


```markdown

## Using QIIME2 Plugin for GCN Normalization

To normalize your data by gene copy number (GCN) using the QIIME2 plugin, follow the steps below.
For more information, visit the q2-gcn-norm GitHub repository (https://github.com/Jiung-Wen/q2-gcn-norm).

### Command
Run the following command to perform GCN normalization:

qiime gcn-norm copy-num-normalize \
  --i-table table-dada2.qza \
  --i-taxonomy taxonomy.qza \
  --o-gcn-norm-table table-normalized.qza
```


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


```r
# Make a new directory and set it as your working directory
create_directory("DspikeIn_16S_OTU", set_working_dir = TRUE)
getwd()

# Please note that these functions have been primarily written based on the 
# phyloseq(https://github.com/joey711/phyloseq) and tidyverse/dplyr (https://github.com/tidyverse/dplyr) packages.
# Therefore, please start by creating a phyloseq object and follow the instructions.
# To create your phyloseq object, please refer to the phyloseq tutorial (https://joey711.github.io/phyloseq).
# The phyloseq object needs to include OTU/ASV, Taxa, phylogenetic tree, DNA reference, 
# and metadata containing spiked species volume, starting from 0 (no spike species added) to 4 (4 μl of spike cell added).

# Briefly:
otu <- read.csv("otu.csv", header = TRUE, sep = ",", row.names = 1)
tax <- read.csv("tax.csv", header = TRUE, sep = ",", row.names = 1)
meta <- read.csv("metadata.csv", header = TRUE, sep = ",")

# Convert data to appropriate formats
meta <- as.data.frame(meta)
taxmat <- as.matrix(tax)
otumat <- as.matrix(otu)
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxmat)

row.names(meta) <- sample_names(OTU)
metadata <- sample_data(meta)
physeq <- phyloseq(OTU, TAX, metadata)
MyTree <- read.tree("tree.nwk")
reference_seqs <- readDNAStringSet(file = "dna-sequences.fasta", format = "fasta")
physeq_16S <- merge_phyloseq(physeq, reference_seqs, MyTree)
physeq_16S <- subset_taxa(physeq_16S, apply(tax_table(physeq_16S), 1, function(x) all(x != "" & !is.na(x))))
physeq_16S <- tidy_phyloseq(physeq_16S)

saveRDS(physeq_16S, file = "physeq_16S.rds")
physeq_16S <- readRDS("physeq_16S.rds")

# Ensure your metadata contains spiked volumes:
# physeq_ITS@sam_data$spiked_volume


```


# Prepare the required information 


```r
# Required Information

# 16S rRNA
spiked_cells <- 1847
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
merged_spiked_species <- "Tetragenococcus_halophilus"
Tetragenococcus_halophilus <- subset_taxa(physeq_16S, Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp")
hashcodes <- row.names(tax_table(Tetragenococcus_halophilus))

# If you intend to use the hashcodes to identify your spiked species, please skip this line of code:
# taxa_names(physeq_16S) <- paste0("ASV", seq(ntaxa(physeq_16S)))

# ITS rDNA
spiked_cells <- 733
species_names <- spiked_species <- merged_spiked_species <- "Dekkera_bruxellensis"
```


# Plot phylogenetic tree with Bootstrap Values
This step will be helpful for handling ASVs with/without Gene Copy Number Correction
This section demonstrates how to use various functions from the package to plot and analyze phylogenetic trees.



```r
# In case there are still several ASVs rooted from the spiked species, you may want to check the phylogenetic distances.
# We first reads DNA sequences from a FASTA file, to perform multiple sequence alignment and compute a distance matrix using the maximum likelihood method, then we construct a phylogenetic tree
# using the Neighbor-Joining method  based on a Jukes-Cantor distance matrix, and plots the tree with bootstrap values.
# we compare the Sanger read of Tetragenococcus halophilus with the FASTA sequence of Tetragenococcus halophilus from our phyloseq object.

# Subset the phyloseq object to include only Tetragenococcus species
Tetragenococcus <- subset_taxa(physeq_16S_ASVs, Species == "Tetragenococcus_halophilus")
Tetragenococcus <- subset_taxa(Tetragenococcus, !is.na(taxa_names(Tetragenococcus)) & taxa_names(Tetragenococcus) != "")
tree <- phy_tree(Tetragenococcus)

# Extract DNA sequences from the phyloseq object
ref_sequences <- refseq(Tetragenococcus)

# Write the sequences to a FASTA format and add the sanger fasta of Tetragenococcus positve control 
writeXStringSet(ref_sequences, "tetra.fasta")

# Now we can run these functions together
# Plot phylogenetic tree
plot_tree(Tetragenococcus, output_prefix = "p0", width = 18, height = 18)

# Plot the tree with glommed OTUs at 0.2 resolution/ or modify it
plot_glommed_tree(Tetragenococcus, resolution = 0.2, output_prefix = "top", width = 18, height = 18)

# Plot the phylogenetic tree with multiple sequence alignment
plot_tree_with_alignment(Tetragenococcus, output_prefix = "tree_alignment", width = 15, height = 15)

# Plot phylogenetic tree with bootstrap values and cophenetic distances
Bootstrap_phy_tree_with_cophenetic(Tetragenococcus, output_file = "tree_with_bootstrap_and_cophenetic.png", bootstrap_replicates = 500)

# Plot Neighbor-Joining tree with bootstrap values
plot_tree_nj("tetra.fasta", output_file = "neighbor_joining_tree_with_bootstrap.png")

```


Figure 1.

| Neighbor Joining Tree with Bootstrap | Tetra Plot with Bootstrap | Cophenetic tree with Bootstrap|
|:------------------------------------:|:----------:|:---------:|
| ![Neighbor Joining Tree](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/JU.png?raw=true) | ![Tetra Plot](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/tetra.png?raw=true) | ![cophenetic](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/Rplot02.png?raw=true) |





```markdown
## Aligned Sequences

The result of the aligned sequences is shown below:
DNAStringSet object of length 5:
    width seq                                                                                                                                         names               
[1]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGG-------------------- 2ddb215ff668b6a24...
[2]   292 GTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGCAAACAGGATTAGATACCCTGGTAGTCC Tetragenococcus h...
[3]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGGCTCGAAAGCGTAGGTAGCAAACAGG-------------------- 65ab824f29da71010...
[4]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGTCTGTAACTGACGCTGAGACTCGAAAGCGTGGGTAGCAAACAGG-------------------- e49935179f23c00fb...
[5]   292 -------------------TACGTAGGTGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGC...CTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAGG-------------------- 0350f990080b4757a...

```



# By now, we have an idea of what works well for our approach, whether it's OTUs or ASVs. 
We will proceed with OTUs. To learn more about why we chose OTUs over ASVs, and what percentage of retrieved spiked species can be considered passed or failed and why, please read our soon-to-be-published paper. For now, check the schematic of our experiment below;

_The VSEARCH with de novo robust clustering algorithms at a 97% similarity threshold was used to reduce potential mistakes._ [Westcott and Schloss 2015](https://doi.org/10.7717/peerj.1487)



![WhyOTUs](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/WhyOTUsOverASVs.png)



## Subsetting and Preprocessing Spiked Data


*Subset the part of the data which is spiked. Keep solely spiked samples using the `spiked_volume` column.*


```r

# Subset spiked samples (264 samples are spiked)
spiked_16S_OTU <- subset_samples(physeq_16S_OTU, spiked_volume %in% c("2", "1"))
spiked_16S_OTU <- tidy_phyloseq(spiked_16S_OTU)

```

### Examine Your Count Data/Biome Before Going Further


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
  threshold_percentage = 0.001,
  threshold_mean_abundance = 0.001,
  threshold_count = 5,
  threshold_relative_abundance = 0.001)
summ_count_phyloseq(FTspiked_16S)

# Random subsampling to even depth with a small trim
spiked_16S_evenDepth <- randomsubsample_Trimmed_evenDepth(spiked_16S_OTU, smalltrim = 0.001)
summ_count_phyloseq(spiked_16S_evenDepth)

```


## Preprocessing for Scaling Factor Calculation  

If you are using OTUs and have only one OTU rooted from the spiked species, you can skip this preprocessing step. Follow the steps below to estimate the success of spike-in, particularly check if you have any samples with under or over-spikes.


**Modify the threshold of acceptable spiked species % as needed. For detailed guidance on acceptable thresholds (passed_range), please refer to the instructions in our upcoming paper.**

```r
# Merg the spiked species
species_name <- "Tetragenococcus_halophilus"
Spiked_16S_OTU_scaled <- Pre_processing_species(physeqASV16, species_name, merge_method = "sum", output_file = "merged_physeq_sum.rds")
Spiked_16S_OTU_scaled <- Pre_processing_species(physeqASV16, species_name, merge_method = "max", output_file = "merged_physeq_max.rds")

Spiked_16S_OTU_scaled <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "sum", output_prefix = "merged_physeq_sum")
Spiked_16S_OTU_scaled <- Pre_processing_hashcodes(physeqASV16, hashcodes, merge_method = "max", output_prefix = "merged_physeq_max")
summ_count_phyloseq(Spiked_16S_OTU_scaled)


Spiked_16S_OTU_scaled <- tidy_phyloseq(Spiked_16S_OTU_scaled)

Spiked_16S_OTU_scaled <- calculate_spike_percentage_hashcodes(Spiked_16S_OTU_scaled, hashcodes, output_path = NULL, passed_range = c(0.1, 10))
Spiked_16S_OTU_scaled <- calculate_spike_percentage_hashcodes(Spiked_16S_OTU_scaled, hashcodes, output_path = NULL, passed_range = c(0.1, 35))
calculate_summary_stats_table(Spiked_16S_OTU_scaled)

merg<-calculate_spike_percentage_species(Spiked_16S_OTU_scaled, spiked_species, identifier_type = "species", output_path = NULL, passed_range = c(0.1, 10))
merg<-calculate_spike_percentage_species(Spiked_16S_OTU_scaled, spiked_species, identifier_type = "species", output_path = NULL, passed_range = c(0.1, 35))
calculate_summary_stats_table(merg)

```

## Data Normalization and Transformation
*Experiment Repetition*

Getting help from [Yerk et al., 2024](https://doi.org/10.1186/s40168-023-01747-z), we checked if we needed to normalize our data before/after calculating our spiked species to account for spiked volume variations and library preparation. We evaluated the need for compositionally aware data transformations, including centered log-ratio (CLR) transformation, DESeq2 variance stabilizing transformation (`run_vst_analysis`), subsampling with a reduced factor for count data (`random_subsample_WithReductionFactor`), proportion adjustment (`proportion.adj`), and prevalence adjustment (`adjusted_prevalence`). Additionally, we considered compositionally naïve data transformations, such as raw data and relative abundance-based transformations (`relativized_filtered_taxa`), before calculating spike-in scaling factors. The only significant variation in the percentage of retrieved spiked species was relevant to VST, so we continued with raw data.


You can repeat the experiment by transforming the data, calculating spike percentage using `calculate_spike_percentage_species()` or `calculate_spike_percentage_hashcodes()`, then checking for the homogeneity of variances using `Bartlett.test()` and ensuring the data is normally distributed using `Shapiro_Wilk_test()`. Finally, plot the results using `transform_plot()`.

```r

methods<-read.csv("methods.csv")

# Check homogeneity of variances
Bartlett.test(methods)

# Check if data is normally distributed
Shapiro_Wilk_test(methods)

y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")

# scale data
scaled <- methods %>% mutate_at(c("Total.reads", "Spike.reads", "Spike.percentage" ), ~(scale(.) %>% as.vector))
 
# Perform Kruskal-Wallis test
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", colors = MG, stat_test = "kruskal.test")

# Perform one-way ANOVA
transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", colors = MG, stat_test = "anova")



```


| Spiked sp Percentage ANOVA | Spiked sp Reads ANOVA | Total Reads ANOVA |
|:--------------------------:|:---------------------:|:-----------------:|
| ![Spike Percentage ANOVA](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/plot_Spike.percentage_ANOVA.jpg?raw=true) | ![Spike Reads ANOVA](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/plot_Spike.reads_ANOVA.png?raw=true) | ![Total Reads ANOVA](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/plot_Total.reads_ANOVA.jpg?raw=true) |

---



### Estimating Scaling Factors After Pre-Processing

To estimate scaling factors, ensure you have the `merged_spiked_species` data, which contains the merged species derived from the spiking process.

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


# Let's check the final conclusion and get the report table of spiked species success or failure.



```r

# Define the parameters once. Do not repeat if already done.
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
identifier_type <- "species"
output_path <- "spike_success_report.docx"

physeq_16S_adj_scaled_perc <- subset_samples(physeq_ITS_adj_scaled_absolute_abundance, sample_or_blank != "blank")
conclusion(physeq_16S_adj_scaled_perc, spiked_species, identifier_type, output_path)

```
*Save your file for later. Stay tuned for the rest: transformation and several visualization methods and displaying the importance of converting relative to absolute abundance in the context of microbial ecology.*


```r

taxa_names(physeq_16S_adj_scaled_absolute_abundance) <- paste0("ASV", seq(ntaxa(physeq_16S_adj_scaled_absolute_abundance)))
physeq_16S_adj_scaled_absolute_abundance <- tidy_phyloseq(physeq_16S_adj_scaled_absolute_abundance)
saveRDS(physeq_16S_adj_scaled_absolute_abundance, "physeq_16S_adj_scaled_absolute_abundance.rds")


# simple barplot of taxonomy abundance
plot_relative <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", group = "Env_broad_scale.x", top = 10, x_axes_font = 10, y_axes_font = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = TRUE, output_prefix = "relativized_abundance_plot")
print(plot_relative)

plot_absolute <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", group = "Env_broad_scale.x", top = 10, x_axes_font = 10, y_axes_font = 10, legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE, output_prefix = "non_relativized_abundance_plot")
print(plot_absolute)

# Check abundance distribution via Ridge Plots before and after converting to absolute abundance

ridgeP_before <- ridge_plot_it(spiked_16S, taxrank = "Family", top_n = 10)
ridgeP_after <- ridge_plot_it(physeq_16S_adj_scaled_absolute_abundance, taxrank = "Family", top_n = 10)


# core_microbiome
plot_core_microbiome_custom(physeq_16S_adj_scaled_absolute_abundance,taxrank = "Phylum")+my_custom_theme()
# core.microbiome is automatically saved in your working directory
core.microbiome <- readRDS("core.microbiome.rds")
taxa_barplot(core.microbiome, target_rank = "Genus", normalize = TRUE, treatment_variable = "animal_type")


# shift to dataframe and plot the abundance of taxa across the factors
meli<-psmelt(physeq_16S_adj_scaled_absolute_abundance)
library(dplyr)
result <- meli %>% 
  dplyr::filter(Abundance > 500) %>% 
  dplyr::select(Abundance, factor1, factor2, factor3, Phylum, Genus)
is_alluvia_form(as.data.frame(meli), axes = 1:6, silent = TRUE)
alluvial_plot <- alluvial_plot(data = meli,axes = c(Abundance, factor1, factor2, factor3, Phylum, Genus),abundance_threshold = 500,silent=TRUE)


# selecting the most important ASVs/OTUs through RandomForest classification
rf_physeq <- RandomForest_selected_ASVs(physeq_16S_adj_scaled_absolute_abundance, response_var = "host.species")
plot_asvs_abundance(rf_physeq, response_var = "host.species", x_var = "ecoregion.III", rank_var = "Phylum")


```
