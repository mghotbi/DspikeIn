---

# DspikeIn

The importance of converting relative to absolute abundance in the context of microbial ecology: Introducing the user-friendly DspikeIn R package.

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

---

## Using QIIME2 Plugin for GCN Normalization

To normalize your data by gene copy number (GCN) using the QIIME2 plugin, follow the steps below. For more information, visit the q2-gcn-norm GitHub repository [here](https://github.com/Jiung-Wen/q2-gcn-norm).

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

## OTUs vs. ASVs Debate
Similar to the debate on the application of OTUs vs. ASVs, which has many contrasting views on their benefits and drawbacks ([Callahan et al., 2017](https://doi.org/10.1038/ismej.2017.119); [Schloss et al., 2021](https://doi.org/10.1128/msphere.00191-21)), there are positive ([Schirrmeister et al., 2012](https://doi.org/10.1186/1471-2180-12-177); [Stoddard et al., 2015](https://doi.org/10.1093/nar/gku1201)) and opposing ([Louca et al., 2018](https://doi.org/10.1038/s43705-023-00266-0); [Gao and Wu et al., 2023](https://doi.org/10.1038/s43705-023-00266-0)) opinions on gene copy number correction for the 16S rRNA marker. Meanwhile, several novelties and modifications have been added to copy number correction to improve accuracy ([Perisin et al., 2016](https://doi.org/10.1038/ismej.2015.161); [Gao and Wu, 2023](https://doi.org/10.1038/s43705-023-00266-0)).

---

## Example Code for Using DspikeIn

Here is an example workflow using the **DspikeIn** package:

### Creating a Phyloseq Object
```r
# Make a new directory and set it as your working directory
create_directory("DspikeIn_16S_OTU", set_working_dir = TRUE)
getwd()

# Read input data
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
```

### Preprocessing and Analysis
```r
# Prepare the required information
spiked_cells <- 1847
species_name <- spiked_species <- c("Tetragenococcus_halophilus", "Tetragenococcus_sp")
merged_spiked_species <- "Tetragenococcus_halophilus"
Tetragenococcus_halophilus <- subset_taxa(physeq_16S, Species == "Tetragenococcus_halophilus" | Species == "Tetragenococcus_sp")
hashcodes <- row.names(tax_table(Tetragenococcus_halophilus))

# Subset the part of the data which is spiked
spiked_16S_OTU <- subset_samples(physeq_16S_OTU, spiked_volume %in% c("2", "1"))
spiked_16S_OTU <- tidy_phyloseq(spiked_16S_OTU)

# Data transformation
readAdj16S <- adjust_abundance_one_third(spiked_16S_OTU, factor = 3)
summ_count_phyloseq(readAdj16S)

# Preprocessing for scaling factor calculation
Spiked_16S_OTU_scaled <- Pre_processing_species(spiked_16S_OTU, species_name)
Spiked_16S_OTU_scaled <- tidy_phyloseq(Spiked_16S_OTU_scaled)

# Calculate spike-in factors
result <- calculate_spikeIn_factors(Spiked_16S_OTU_scaled, spiked_cells, merged_spiked_species)
scaling_factors <- result$scaling_factors

# Convert relative counts to absolute counts
physeq_16S_adj_scaled_AbsoluteCount <- convert_to_absolute_counts(Spiked_16S_OTU_scaled, scaling_factors)
```

### Visualizing Data
```r
# Plot relative abundance
plot_relative <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", group = "Env_broad_scale.x", top = 10, output_prefix = "relativized_abundance_plot")
print(plot_relative)

# Plot absolute abundance
plot_absolute <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", group = "Env_broad_scale.x", top = 10, output_prefix = "non_relativized_abundance_plot")
print(plot_absolute)
```

---

## Conclusion

By now, we have an idea of what works well for our approach, whether it's OTUs or ASVs. To learn more about why we chose OTUs over ASVs, and what percentage of retrieved spiked species can be considered passed or failed and why, please read our soon-to-be-published paper. 

For more details and advanced usage, stay tuned for our publication and additional documentation.

---

Contact: Mitra Ghotbi [mitra.ghotbi@gmail.com](mailto:mitra.ghotbi@gmail.com)

GitHub Repository: [DspikeIn](https://github.com/mghotbi/DspikeIn.git)

---