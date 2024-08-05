install_and_load <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("phyloseq", "DESeq2", "edgeR", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics", "ShortRead", "Biostrings", "DECIPHER", "microbiome")) {
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

install_required_packages <- function() {
  cran_packages <- c(
    "dplyr", "ape", "ggplot2", "ggalluvial", "phangorn", "data.table", "matrixStats", "ggridges", "flextable", "RColorBrewer", "randomForest", "ggtree", "ggpubr", "agricolae"
  )
  
  bioc_packages <- c(
    "phyloseq", "microbiome", "DESeq2", "ANCOMBC", "ShortRead", "msa", "Biostrings", "DECIPHER", "edgeR", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan"
  )
  
  github_packages <- c(
    "mikemc/speedyseq", "microsud/microbiomeutilities", "vmikk/metagMisc"
  )
  
  installed <- installed.packages()
  
  to_install_cran <- cran_packages[!cran_packages %in% installed[, "Package"]]
  to_install_bioc <- bioc_packages[!bioc_packages %in% installed[, "Package"]]
  
  if (length(to_install_cran) > 0) {
    tryCatch({
      install.packages(to_install_cran, dependencies = TRUE)
    }, error = function(e) {
      message("Error installing CRAN packages: ", e$message)
    })
  }
  
  if (length(to_install_bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    tryCatch({
      BiocManager::install(to_install_bioc, dependencies = TRUE)
    }, error = function(e) {
      message("Error installing Bioconductor packages: ", e$message)
    })
  }
  
  if (length(github_packages) > 0) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    for (pkg in github_packages) {
      tryCatch({
        devtools::install_github(pkg)
      }, error = function(e) {
        message("Error installing GitHub package ", pkg, ": ", e$message)
      })
    }
  }
  
  # Load all packages
  lapply(c(cran_packages, bioc_packages, gsub(".*/", "", github_packages)), function(pkg) {
    tryCatch({
      library(pkg, character.only = TRUE)
    }, error = function(e) {
      message("Error loading package ", pkg, ": ", e$message)
    })
  })
  
  message("All required packages are installed and loaded.")
}

# Run the function to install and load all packages
install_required_packages()
