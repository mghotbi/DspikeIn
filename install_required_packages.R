# Function to install and load required packages
install_and_load <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("phyloseq", "DESeq2", "edgeR", "PoissonSeq", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics", "ShortRead", "Biostrings", "DECIPHER", "microbiome")) {
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

# Specific installation for PoissonSeq from source if needed
if (!requireNamespace("PoissonSeq", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/PoissonSeq_1.1.2.tar.gz", repos = NULL, type = "source")
}

# List of required packages
required_packages <- c(
  "phyloseq", "DESeq2", "edgeR", "PoissonSeq", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics", "ShortRead", "Biostrings", "DECIPHER", "microbiome",
  "dplyr", "ape", "ggplot2", "ggalluvial", "phangorn", "msa", "data.table", "matrixStats", "ggridges", "flextable", "RColorBrewer", "randomForest", "ggtree", "ggpubr"
)

# Install and load required packages for DspikeIn
install_required_packages <- function() {
  cran_packages <- c(
    "dplyr", "ape", "ggplot2", "ggalluvial", "phangorn", "msa", "data.table", "matrixStats", "ggridges", "flextable", "RColorBrewer", "randomForest", "ggtree", "ggpubr"
  )
  
  bioc_packages <- c(
    "phyloseq", "microbiome", "DESeq2", "ShortRead", "Biostrings", "DECIPHER", "edgeR", "PoissonSeq", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics"
  )
  
  github_packages <- c(
    "mikemc/speedyseq", "microsud/microbiomeutilities", "vmikk/metagMisc"
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

# Run the function to install and load all packages
install_required_packages()
