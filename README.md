# DspikeIn
The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package

![How it works](https://github.com/mghotbi/DspikeIn/blob/MitraGhotbi/DspikeIn.png)


---

## Introduction

In our study, *Tetragenococcus halophilus* and *Dekkera bruxellensis* were selected as taxa to spike into gut microbiome samples  based on our previous studies.
## Methodology

### 1. Growth of Stock Cell Suspensions
- **Tetragenococcus halophilus**: Cultivated in tryptic soy broth.
- **Dekkera bruxellensis**: Cultivated in potato dextrose broth.
- Both microbial cultures were serially diluted, and optical density (OD) measurements were obtained using a ClarioStar plate reader.

### 2. DNA Extraction

- DNA was extracted using the Qiagen DNeasy Powersoil Pro Kit.
- These DNA isolations served as standards to determine the appropriate spike-in volume of cells to represent 0.1-10% of a sample, as detailed in [this Nature article](https://www.nature.com/articles/s41586-021-03241-8).

---

## Bioinformatics Processing and Downstream Analyses

### Gene Marker Analysis
- 16S rRNA and ITS rDNA gene markers were analyzed.

### Normalization
- Normalization was performed for community-weighted mean ribosomal operon copy numbers.

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


Similar to the debate on the application of OTUs vs. ASVs, which has many contrasting views on their benefits and drawbacks ([Callahan et al., 2017](https://doi.org/10.1038/ismej.2017.119); [Schloss et al., 2021](https://doi.org/10.1128/msphere.00191-21)), there are differing opinions on estimating gene copy numbers and making corrections accordingly. Positive views ([Gao and Wu et al., 2023](https://doi.org/10.1038/s43705-023-00266-0)) contrast with negative perspectives ([Louca et al., 2018](https://doi.org/10.1038/s43705-023-00266-0)).


Using our spike-in positive controls and assessing the percentage of retrieved species, we were able to compare the results at each step before selecting a pathway. We utilized the method by [Gao and Wu et al. (2023)](https://github.com/wu-lab-uva/16S-rRNA-GCN-Predcition) to normalize for 16S rRNA gene marker copy numbers. Due to the variability in rDNA gene copy numbers, straightforward translation of rDNA read counts into the abundance of individual organisms is precluded ([Lavrinienko et al., 2021](https://doi.org/10.1016/j.tim.2020.05.019)).


For ITS, we did not need to use copy number correction. However, we recommend conducting a literature review before deciding whether to sum or select the maximum abundance of OTUs/ASVs rooted from spiked-in species and calculating spike-in factors. We believe systematic evaluation before selecting or refuting each method can help prevent miscalculations ([Lofgren et al., 2018](https://doi.org/10.1111/mec.14995)), which can aid in establishing a system-dependent method for copy number correction in ITS markers.

---
```r
# Make a new directory and set it as your working directory
create_directory("ExampleITS", set_working_dir = TRUE)
getwd()

# Please note that these functions have been primarily written based on the 
# phyloseq(https://github.com/joey711/phyloseq) and microbiome (https://github.com/microbiome/microbiome) packages.
# Therefore, please start by creating a phyloseq object and follow the instructions.
# To create your phyloseq object, please refer to the phyloseq tutorial (https://joey711.github.io/phyloseq/import-data.html).
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

# Alternatively, you can make your phyloseq object using the BIOM JSON format:
# mydata <- import_biom(BIOMfilename = "taxonomy.feature-table.biom", parseFunction = parse_taxonomy_default)
# mapfile <- import_qiime_sample_data("metadata.txt")
# physeq_16S <- merge_phyloseq(mydata, mapfile, tree, reference_seqs)
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
