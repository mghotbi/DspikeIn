#' Install and Load Necessary Packages
#'
#' This function checks if the required packages are installed, installs any missing packages,
#' and loads them into the R session.
#'
#' @return NULL. This function is called for its side effects of installing and loading packages.
#' @examples
#' # Install and load necessary packages
#' install_load_packages_helper()
#' @export
install_load_packages_helper <- function() {
  reqpkg <- c(
    "edgeR", "microbiome", "matrixStats", "e1071", "vcd", "qpdf",
    "optparse", "phyloseq", "caret", "DESeq2", "foreach", "magrittr",
    "doParallel", "viridis", "plyr", "reshape2", "ggplot2", "grid",
    "vegan", "scales", "cluster", "ape", "dplyr", "ggridges", "speedyseq",
    "microbiomeutilities", "MicEco", "intergraph", "ggalluvial", "phangorn",
    "msa", "data.table", "decontam", "ggtext", "devtools", "dada2",
    "ggtree", "Biostrings", "ggpubr", "agridat", "rstatix", "emmeans",
    "randomForest", "DECIPHER", "reshape", "agricolae", "gridExtra",
    "flextable", "RColorBrewer", "microbial", "descr", "ShortRead"
  )
  
  for (pkg in reqpkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}
