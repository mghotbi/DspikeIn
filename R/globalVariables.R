#' @title title Suppress R CMD check warnings for global variables
#'
#' @description This file uses `utils::globalVariables()` to declare global variables in order to suppress warnings
#' from R CMD check about undefined global variables used in non-standard evaluation (NSE).
#'
#' @name globalVariables_suppress
#' @return NULL. Used to suppress global variable warnings in R CMD check.
#' @noRd
#' @importFrom utils globalVariables write.csv install.packages capture.output
#' @importFrom ggstar geom_star
#' @importFrom grDevices rainbow
#' @importFrom phyloseq phyloseq taxa_sums tax_table sample_sums sample_data tax_glom prune_taxa subset_taxa transform_sample_counts taxa_names psmelt taxa_names<- otu_table<-
#' @importFrom ggplot2 ggplot aes aes_string geom_bar geom_label geom_point geom_text ggtitle ggsave scale_fill_manual scale_color_manual theme theme_minimal element_blank element_line element_rect element_text guide_legend xlab ylab
#' @importFrom dplyr filter mutate pull summarise group_by ungroup desc all_of top_n select rename bind_cols
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix estimateSizeFactorsForMatrix results
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom randomForest randomForest importance
#' @importFrom ape boot.phylo
#' @importFrom ggtree ggtree theme_tree2 geom_tiplab geom_tippoint geom_text2
#' @importFrom stats wilcox.test median quantile sd p.adjust reformulate setNames
#' @importFrom flextable save_as_docx flextable fontsize font color bold italic save_as_docx
#' @importFrom Biostrings readDNAStringSet
#' @importFrom phangorn phyDat dist.ml
#' @importFrom msa msa msaConvert
#' @importFrom microbiome plot_core
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics mtext boxplot
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 after_stat ggplot aes geom_label theme scale_x_discrete scale_fill_manual ylab ggtitle guides guide_legend facet_grid
#' @importFrom ggpubr color_palette
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom DECIPHER AlignSeqs
#' @importFrom scales label_scientific
#' @importFrom ggalluvial is_alluvia_form geom_alluvium geom_stratum
#' @importFrom grid unit
#' @importFrom limma makeContrasts
#' @importFrom ggpubr stat_regline_equation stat_cor
#' @importFrom SummarizedExperiment colData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment rowTree
#' @importFrom S4Vectors metadata
#' @importFrom stats aggregate reorder
#' @importFrom dplyr where
#' @importFrom data.table :=
#' @importFrom methods is
utils::globalVariables(c(
  ".data", "Taxon_Label", "prevalence", "Metrics", "Value", "Label", "log_mean_abundance", "LR", "Genus", "lfcSE", "Degree", "prune_common_taxa", "reorder", "LFC_Direction", "read_spike_report", "Among_Module_Connectivity", "EigenvectorCentrality",
  "Within_Module_Connectivity", "group", "eq.label", "tail", "PageRank", "Closeness", "Betweenness", "Community",
  "Strength", "Local_Efficiency", "Coreness", "Redundancy", "get_refseq",
  "Constraint", "Efficiency", "Z_Score", "Significance", "format_common_data", "FDR", "contrast", "model.matrix", "convert_phyloseq_to_tse", "aggregateByTaxonomy", "agglomerateByRank", "group_var",
  "taxa_sums", "perform_edgeR", "perform_DESeq2", "diff_abn", "glom_taxa_at_rank", "Class", "..eq.label..", "..rr.label..", "..p.label..", "Sample", "Family", "Species",
  "Abundance", "otu_table<-", "taxa_names<-", "visualize_differential_abundance", "model.matrix", "convert_phyloseq_to_tse", "aggregateByTaxonomy",
  "final_results", "group_var", "point_size", "palette", "StatStratum", "makeContrasts", "treatment", "padj", "agglomerateByRank",
  "DESeq", "DESeqDataSetFromMatrix", "DGEList", "FDR", "MeanDecreaseGini", "Total_Reads_total", "pseudocount", "metadata", "phy_tree",
  "OTU", "PValue", "Percentage", "Result", "TotalAbundance", "Total_Reads_spiked", "assays", "rowData", "colData", "rowTree",
  "aes", "aes_string", "annotate", "as.formula", "boot.phylo", "boxplot", "capture.output", ".", "format_common_data",
  "complete.cases", "counts", "desc", "dev.off", "dist.ml", "element_blank", "ASV", "x", "y", "ecount", "read.csv", "get_reference_seq",
  "element_line", "element_rect", "element_text", "estimateDisp", "estimateGLMCommonDisp", "aggregate", "where",
  "estimateGLMTagwiseDisp", "estimateSizeFactors", "estimateSizeFactorsForMatrix", "RelativeAbundance", "Spiked_Reads", "Total_Reads", "Quadrant",
  "Node", "Step", "LargestComponentRatio", "MG", "filter", "filter_taxa", "geom_bar", "geom_label", "geom_point", "geom_text2",
  "geom_tiplab", "geom_tippoint", "ggplot", "ggsave", "ggtitle", "ggtree", "glmFit", "aggregate",
  "glmLRT", "group_by", "guide_legend", "head", "importance", "install.packages", "as_tibble",
  "isTip", "kruskal", "labs", "log2FoldChange", "logFC", "median", "merge_taxa", ":=", "group_label", "geom_tile",
  "model.matrix", "mutate", "na.omit", "nodelabels", "normCounts", "nsamples", "ntaxa", "%>%",
  "one_of", "otu_table", "par", "phyDat", "phy_tree", "phyloseq", "png", "prune_samples",
  "prune_taxa", "pull", "quantile", "randomForest", "rarefy_even_depth", "readDNAStringSet", "recordPlot",
  "refseq", "resid", "residuals", "results", "rowMedians", "sample_data", "sample_names",
  "sample_sums", "save_as_docx", "scale_color_manual", "scale_fill_manual", "scale_x_discrete",
  "scale_y_continuous", "sd", "stat_compare_means", "str", "subset_taxa", "summarise",
  "tax_glom", "tax_table", "taxa_names", "taxa_sums", "taxon", "theme", "theme_minimal",
  "theme_tree2", "topTags", "transform_sample_counts", "ungroup", "unit", "vars", "wilcox.test",
  "write.csv", "xlab", "ylab", ".data", "taxa_sums", "Class", "Family", "Species", "Abundance",
  "pvalue", "branch.length", "stratum", "label", "node", "visualize_differential_abundance"
))
