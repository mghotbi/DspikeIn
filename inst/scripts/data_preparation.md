# Data Preparation and Provenance for `inst/extdata/` and `data`

This document outlines the creation of synthetic and original test datasets included in the **DspikeIn** package. It describes how `phyloseq` and `TreeSummarizedExperiment` (TSE) objects were derived from raw tabular data and includes examples for generating co-occurrence networks. These datasets support demonstration, testing, and reproducibility of workflows for absolute microbial abundance estimation.

---
  
  ## Dataset Description
  inst/extdata
- **`physeq`** — synthetic `phyloseq` object including spike-in taxa for use in controlled community workflows.  
- **`tse`** — `TreeSummarizedExperiment` equivalent of the synthetic dataset.  
- **`Sample.fasta`** — short read of *Tetragenococcus halophilus* from 16S rRNA sequencing in a herptile microbiome study.
  data
- **`physeq_16SOTU`** — `phyloseq` object built from herptile gut microbiome 16S rRNA sequencing data.  
- **`physeq_ITSOTU`** — `phyloseq` object built from herptile gut microbiome ITS1 region sequencing data.

All datasets were generated using reproducible R-based workflows and include taxonomic annotations, OTU tables, and sample metadata.

---
  
## Source and Archiving

- **Author**: Mitra Ghotbi  
- **Institution**: [Herptile Microbiomes](https://herptilemicrobiomes.org/research/)  
- **Funding**: National Science Foundation (NSF) grants EF-2125065, EF-2125066, EF-2125067  
- **Manuscript**: Under review at *ISME Journal*  
- **Preprint**: Available at [bioRxiv](https://doi.org/10.1101/2024.12.27.630554)

For full methodological and biological context, see:

Ghotbi, M. *et al.* (2024). *DspikeIn: A spike-in controlled absolute abundance framework for microbial community analysis*. *bioRxiv*. https://doi.org/10.1101/2024.12.27.630554


---
  
  ## SRA Accessions
  
  Raw sequencing data are publicly available at the NCBI Sequence Read Archive:
  
- [PRJNA1202922](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1202922)  
- [PRJNA1210664](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1210664)

---
  
  ## Preprocessing Overview
  
  All 16S rRNA and ITS1 datasets were processed using **QIIME2 v2024.10** as follows:
  
  - **Sequencing Platform**: Illumina MiSeq  
- **Demultiplexing & Primer Removal**: Cutadapt  
- **Paired-End Merging**: FLASH v1.2.11 (truncated to 225 bp forward and 220 bp reverse)  
- **OTU Clustering**: De novo at 97% similarity  
- **Taxonomy Assignment**:  
  - 16S: sklearn classifier with the SILVA 138 database (515F/806R region)  
- ITS1: Region trimming with ITSxpress and classification using the UNITE dynamic database (v10)

Reference databases:  
  [https://github.com/mghotbi/QIIME-compatible-database](https://github.com/mghotbi/QIIME-compatible-database)

---
  
  ## Data Generation Code and Downstream Processing
  
  ### Build `physeq` (Synthetic Community)
  
```r
library(phyloseq)

otu <- read.csv("otu.csv", row.names = 1)
tax <- read.csv("tax.csv", row.names = 1)
meta <- read.csv("metadata.csv")

colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
row.names(meta) <- sample_names(OTU)
metadata <- sample_data(meta)

physeq <- phyloseq(OTU, TAX, metadata)
physeq <- DspikeIn::tidy_phyloseq_tse(physeq)
saveRDS(physeq, "physeq.rds")

tse <- DspikeIn::convert_phyloseq_to_tse(physeq)
saveRDS(tse, "tse.rds")


# Build physeq_16SOTU and physeq_ITSOTU (Herptile Microbiome)

# Load data
otu <- read.csv("otu_table_ITS.csv", sep = ",", row.names = 1, header = TRUE)
tax <- read.csv("taxonomy_table.csv", sep = ",", row.names = 1, header = TRUE)
meta <- read.csv("sample_metadata.csv", row.names = 1)
tree <- ape::read.tree("phy_tree.nwk")
repseqs <- Biostrings::readDNAStringSet("rep_seqs.fasta")

# 16S dataset
physeq_16SOTU <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(meta),
  phy_tree(tree),
  refseq(repseqs)
)
physeq_16SOTU <- DspikeIn::tidy_phyloseq_tse(physeq_16SOTU)
saveRDS(physeq_16SOTU, "physeq_16SOTU.rds")

# ITS dataset
otu <- read.csv("otu_table_ITS.csv", sep = ",", row.names = 1, header = TRUE)
tax <- read.csv("taxonomy_table.csv", sep = ",", row.names = 1, header = TRUE)
meta <- read.csv("sample_metadata.csv", row.names = 1)
tree <- ape::read.tree("phy_tree.nwk")
repseqs <- Biostrings::readDNAStringSet("rep_seqs.fasta")

physeq_ITSOTU <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(meta),
  phy_tree(tree),
  refseq(repseqs)
)
physeq_ITSOTU <- DspikeIn::tidy_phyloseq_tse(physeq_ITSOTU)
saveRDS(physeq_ITSOTU, "physeq_ITSOTU.rds")


