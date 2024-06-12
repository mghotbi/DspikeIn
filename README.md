# DspikeIn
The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package

![How it works](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/DspikeIn.png)


---

## Introduction

In our study, *Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples  based on our previous studies.
## Methodology

### Growth of Stock Cell Suspensions
- **Tetragenococcus halophilus**: Cultivated in tryptic soy broth.
- **Dekkera bruxellensis**: Cultivated in potato dextrose broth.
- Both microbial cultures were serially diluted, and optical density (OD) measurements were obtained using a ClarioStar plate reader.

### DNA Extraction

- DNA was extracted using the Qiagen DNeasy Powersoil Pro Kit.
- These DNA isolations served as standards to determine the appropriate spike-in volume of cells to represent 0.1-10% of a sample, as detailed in [Roa et al., 2021](https://www.nature.com/articles/s41586-021-03241-8).

---

## Bioinformatics Processing and Downstream Analyses

### Gene Marker Analysis
- 16S rRNA and ITS rDNA gene markers were analyzed.

### Normalization
- Normalization was performed for community-weighted mean ribosomal operon copy numbers.


```markdown

## Using QIIME2 Plugin for GCN Normalization

To normalize your data by gene copy number (GCN) using the QIIME2 plugin, follow the steps below.
For more information, visit the q2-gcn-norm GitHub repository (https://github.com/Jiung-Wen/q2-gcn-norm).

### Command
Run the following command to perform GCN normalization:

```bash
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



Using our spike-in positive controls and assessing the percentage of retrieved spiked species, we were able to compare the results at each step before selecting a pathway. We utilized q2-gcn-norm based on rrnDB database (version 5.7) to normalize for 16S rRNA gene marker copy numbers ([qiime2 plugin; gcn-norm](https://github.com/Jiung-Wen/q2-gcn-norm)). Due to the variability in rDNA gene copy numbers, straightforward translation of rDNA read counts into the abundance of individual organisms is precluded ([Lavrinienko et al., 2021](https://doi.org/10.1016/j.tim.2020.05.019)).



In fact, for ITS we did not need to use copy number correction. However, we recommend conducting a literature review before deciding whether to sum or select the maximum abundance of OTUs/ASVs rooted from spiked-in species and calculating spike-in factors. We believe systematic evaluation before selecting or refuting each method can help prevent miscalculations ([Lofgren et al., 2018](https://doi.org/10.1111/mec.14995)), which can aid in establishing a system-dependent method for copy number correction in ITS markers.

---
```r
# Make a new directory and set it as your working directory
create_directory("ExampleITS", set_working_dir = TRUE)
getwd()

# Please note that these functions have been primarily written based on the 
# phyloseq(https://github.com/joey711/phyloseq) and microbiome (https://github.com/microbiome/microbiome) packages.
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




```r
# This step is helpful for handling ASVs with/without Gene Copy Number Correction
# This section demonstrates how to use various functions from the package to plot and analyze phylogenetic trees.
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
![neighbor_joining_tree_with_bootstrap](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/trees.png)




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
Moreover, we have learned which ASVs/OTUs need to be merged and whether their abundance should be summed up or if the maximum abundance is sufficient. We will proceed with OTUs. We used VSEARCH with de novo robust clustering algorithms at a 97% similarity threshold to reduce potential mistakes.



```markdown

### Preprocessing for ASVs

If you are using OTUs and have only one OTU rooted from the spiked species, you can skip this pre-processing step. For those using ASVs, follow the steps below:

### Examine Your Count Data/Biome Before Calculating Absolute going further
# Summarize the initial statistics for ASVs
initial_stat_ASV <- summ_phyloseq_ASVID(physeq_16S_OTU)

# Summarize the initial statistics sample-wise
initial_stat_sampleWise <- summ_phyloseq_sampleID(physeq_16S_OTU)

# Summarize the count data
summ_count_phyloseq(physeq_16S_OTU)

# Check the summary statistics
# Ensure the input is in dataframe format for this function
calculate_summary_stats_table(initial_stat_sampleWise)

```
