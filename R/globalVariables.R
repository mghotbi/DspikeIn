#' Suppress R CMD check warnings for global variables
#'
#' This file uses `utils::globalVariables()` to declare global variables in order to suppress warnings
#' from R CMD check about undefined global variables used in non-standard evaluation (NSE).
#'
#' @name globalVariables  # Add the @name tag to name this block
#' @importFrom utils globalVariables
#' @importFrom phyloseq phyloseq taxa_sums tax_table sample_sums sample_data
#' @importFrom ggplot2 ggplot aes aes_string geom_bar geom_label geom_point geom_text ggtitle ggsave scale_fill_manual scale_color_manual theme theme_minimal element_blank element_line element_rect element_text guide_legend xlab ylab
#' @importFrom dplyr filter mutate pull summarise group_by ungroup desc 
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix estimateSizeFactorsForMatrix results
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom randomForest randomForest importance
#' @importFrom ape boot.phylo
#' @importFrom ggtree ggtree theme_tree2 geom_tiplab geom_tippoint geom_text2
#' @importFrom stats wilcox.test median quantile sd
#' @importFrom utils write.csv install.packages capture.output
#' @importFrom flextable save_as_docx
#' @importFrom Biostrings readDNAStringSet
utils::globalVariables(c(
  ".data", "taxa_sums", "Class", "Family", "Species", "Abundance",
  "DESeq", "DESeqDataSetFromMatrix", "DGEList", "FDR", "MeanDecreaseGini", 
  "OTU", "PValue", "Percentage", "Result", "TotalAbundance", "Total_Reads_spiked",
  "aes", "aes_string", "annotate", "as.formula", "boot.phylo", "boxplot", "capture.output", 
  "complete.cases", "counts", "desc", "dev.off", "dist.ml", "element_blank", 
  "element_line", "element_rect", "element_text", "estimateDisp", "estimateGLMCommonDisp",
  "estimateGLMTagwiseDisp", "estimateSizeFactors", "estimateSizeFactorsForMatrix", 
  "filter", "filter_taxa", "geom_bar", "geom_label", "geom_point", "geom_text2", 
  "geom_tiplab", "geom_tippoint", "ggplot", "ggsave", "ggtitle", "ggtree", "glmFit",
  "glmLRT", "group_by", "guide_legend", "head", "importance", "install.packages", 
  "is", "isTip", "kruskal", "labs", "log2FoldChange", "logFC", "median", "merge_taxa", 
  "model.matrix", "mutate", "na.omit", "nodelabels", "normCounts", "nsamples", "ntaxa", 
  "one_of", "otu_table", "par", "phyDat", "phy_tree", "phyloseq", "png", "prune_samples", 
  "prune_taxa", "pull", "quantile", "randomForest", "rarefy_even_depth", "readDNAStringSet", 
  "refseq", "resid", "residuals", "results", "rowMedians", "sample_data", "sample_names", 
  "sample_sums", "save_as_docx", "scale_color_manual", "scale_fill_manual", "scale_x_discrete",
  "scale_y_continuous", "sd", "stat_compare_means", "str", "subset_taxa", "summarise", 
  "tax_glom", "tax_table", "taxa_names", "taxa_sums", "taxon", "theme", "theme_minimal", 
  "theme_tree2", "topTags", "transform_sample_counts", "ungroup", "unit", "vars", "wilcox.test", 
  "write.csv", "xlab", "ylab",".data", "taxa_sums", "Class", "Family", "Species", "Abundance", 
  "pvalue", "branch.length", "stratum", "label", "node"
))
