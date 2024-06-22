#' Install and load required packages for YourPackageName
#'
#' This function installs and loads the packages required by YourPackageName.
#'

install_required_packages <- function() {
  cran_packages <- c(
    "dplyr",
    "ape",
    "ggplot2",
    "ggalluvial",
    "phangorn",
    "msa",
    "data.table",
    "matrixStats",
    "ggridges",
    "agricolae",
    "flextable",
    "RColorBrewer",
    "randomForest",
    "ggtree",
    "ggpubr",
    "methods",
    "tidyverse"
  )
  
  bioc_packages <- c(
    "phyloseq",
    "microbiome",
    "microbiomeutilities",
    "speedyseq",
    "DESeq2",  
    "ShortRead",
    "Biostrings", 
    "DECIPHER"
  )
  
  github_packages <- c(
    "microbiome/microbiome",
    "vmikk/metagMisc"
  )
  
  installed <- installed.packages()
  
  to_install_cran <- cran_packages[!cran_packages %in% installed[, "Package"]]
  to_install_bioc <- bioc_packages[!bioc_packages %in% installed[, "Package"]]
  
  if (length(to_install_cran) > 0) {
    install.packages(to_install_cran, dependencies = TRUE)
  }
  
  if (length(to_install_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(to_install_bioc, dependencies = TRUE)
  }
  
  if (length(github_packages) > 0) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    for (pkg in github_packages) {
      devtools::install_github(pkg)
    }
  }
  
  # Load all packages
  lapply(c(cran_packages, bioc_packages, gsub(".*/", "", github_packages)), library, character.only = TRUE)
  
  message("All required packages are installed and loaded.")
}
