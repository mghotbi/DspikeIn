#' Perform and Visualize Differential Abundance Analysis with edgeR or DESeq2
#'
#' This function performs differential abundance analysis using either edgeR or DESeq2 
#' and generates a volcano plot based on microbiome data.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param method A character string specifying the method to use: "edgeR" or "DESeq2".
#' @param group_var A string specifying the grouping variable in sample data.
#' @param contrast A vector specifying the levels to compare (e.g., c("Control", "Treated")).
#' @param pseudocount A numeric value for pseudocount addition to handle zeros. Default is 1.
#' @param significance_level A numeric value specifying the significance threshold (FDR). Default is 0.05.
#' @param output_csv_path Path to save the results as a CSV file. Default is "DA_results.csv".
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 3.
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa. Default is "Genus".
#' @param palette A character vector of color hex codes for plotting. Default is `color_palette$extended_palette`.
#' @return A list containing the analysis results, a phyloseq object with significant OTUs/ASVs, and a volcano plot.
#'
#' @importFrom phyloseq otu_table sample_data tax_table prune_taxa tax_glom
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma makeContrasts
#' @importFrom dplyr arrange mutate filter pull left_join rename
#' @importFrom ggplot2 ggplot aes geom_point scale_y_continuous theme_minimal labs scale_color_manual scale_shape_manual theme
#' @importFrom scales label_scientific
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' # Using edgeR
#' results_edgeR <- perform_and_visualize_DA(
#'   ps = ps,
#'   method = "edgeR",
#'   group_var = "Treatment",
#'   contrast = c("Control", "Treated"),
#'   output_csv_path = "DA_edgeR.csv",
#'   target_glom = "Genus",
#'   significance_level = 0.05
#' )
#' 
#' # Using DESeq2
#' results_DESeq2 <- perform_and_visualize_DA(
#'   ps = ps,
#'   method = "DESeq2",
#'   group_var = "Treatment",
#'   contrast = c("Control", "Treated"),
#'   output_csv_path = "DA_DESeq2.csv",
#'   target_glom = "Genus",
#'   significance_level = 0.05
#' )
#'
#' # Visualize Volcano Plot
#' print(results_edgeR$plot)
#' ps_sig <-results_edgeR$ps_significant
#' print(results_DESeq2$plot)
#' }
#' 
#' @export
perform_and_visualize_DA <- function(ps, method, group_var, contrast, pseudocount = 1, significance_level = 0.05, 
                                     output_csv_path = "DA_results.csv", point_size = 3, 
                                     target_glom = "Genus", palette = color_palette$extended_palette) {
  if (!method %in% c("edgeR", "DESeq2")) {
    stop("Invalid method. Please choose 'edgeR' or 'DESeq2'.")
  }
  
  # Step 1: Preprocess phyloseq object
  ps <- remove_zero_negative_count_samples(ps, pseudocount = pseudocount)
  ps <- convert_categorical_to_factors(ps)
  ps <- glom_taxa_at_rank(ps, rank = target_glom)
  
  # Step 2: Perform analysis
  if (method == "edgeR") {
    results <- perform_edgeR(ps, group_var, contrast, significance_level)
  } else {
    results <- perform_DESeq2(ps, group_var, contrast, significance_level)
  }
  
  # Step 3: Save results to CSV
  utils::write.csv(results, output_csv_path, row.names = FALSE)
  cat(method, "results saved to:", output_csv_path, "\n")
  
  # Step 4: Filter significant OTUs
  significant_otus <- results %>% dplyr::filter(FDR < significance_level) %>% dplyr::pull(OTU)
  ps_significant <- phyloseq::prune_taxa(significant_otus, ps)
  
  # Step 5: Merge results with taxonomy metadata
  tax_metadata <- as.data.frame(phyloseq::tax_table(ps_significant))
  tax_metadata$OTU <- rownames(tax_metadata)
  merged_results <- dplyr::left_join(results, tax_metadata, by = "OTU")
  
  # Step 6: Generate volcano plot
  merged_results$treatment <- ifelse(merged_results$logFC > 0, contrast[1], contrast[2])
  p <- ggplot2::ggplot(merged_results, ggplot2::aes(x = logFC, y = -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(color = diff_abn, shape = treatment), size = point_size) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Log2 Fold Change", 
      y = "-log10(P-value)", 
      title = paste("Volcano Plot -", group_var)
    ) +
    ggplot2::theme(legend.position = "bottom")
  
  return(list(results = merged_results, ps_significant = ps_significant, plot = p))
}

# Supporting Functions ---------------------------------------------------

#' Perform edgeR Analysis
perform_edgeR <- function(ps, group_var, contrast, significance_level = 0.05) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  metadata[[group_var]] <- as.factor(metadata[[group_var]])
  
  if (!all(contrast %in% levels(metadata[[group_var]]))) {
    stop("One or both levels in the contrast are not in the group_var factor levels.")
  }
  
  dge <- edgeR::DGEList(counts = otu, group = metadata[[group_var]])
  design <- stats::model.matrix(~ 0 + metadata[[group_var]])
  colnames(design) <- make.names(levels(metadata[[group_var]]))
  
  contrast_matrix <- limma::makeContrasts(contrasts = paste0(make.names(contrast[1]), "-", make.names(contrast[2])), 
                                          levels = colnames(design))
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit, contrast = contrast_matrix)
  
  res <- edgeR::topTags(lrt, n = Inf)$table
  res <- res %>%
    dplyr::mutate(FDR = p.adjust(PValue, method = "BH"), diff_abn = FDR < significance_level) %>%
    dplyr::rename(logFC = logFC, pvalue = PValue) %>%
    dplyr::mutate(OTU = rownames(res))
  
  return(res)
}

#' Perform DESeq2 Analysis
perform_DESeq2 <- function(ps, group_var, contrast, significance_level = 0.05) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  metadata[[group_var]] <- as.factor(metadata[[group_var]])
  
  if (!all(contrast %in% levels(metadata[[group_var]]))) {
    stop("One or both levels in the contrast are not in the group_var factor levels.")
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, design = as.formula(paste("~", group_var)))
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c(group_var, contrast[1], contrast[2]))
  res <- as.data.frame(res)
  
  res <- res %>%
    dplyr::rename(logFC = log2FoldChange) %>%
    dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"), diff_abn = FDR < significance_level) %>%
    dplyr::mutate(OTU = rownames(res))
  
  return(res)
}

# # Usage Example
# results_edgeR <- perform_and_visualize_DA(
#   ps = ps, 
#   method = "edgeR", 
#   group_var = "Treatment", 
#   contrast = c("Control", "Flooding"), 
#   output_csv_path = "DA_edgeR.csv", 
#   target_glom = "Genus", 
#   significance_level = 0.05
# )
# 
# print(results_edgeR$plot)
# head(results_edgeR$results)  # View significant taxa
# results_edgeR$ps_significant
# 
# 
# results_DESeq2 <- perform_and_visualize_DA(
#   ps = ps, 
#   method = "DESeq2", 
#   group_var = "Treatment", 
#   contrast = c("Control", "Flooding"), 
#   output_csv_path = "DA_DESeq2.csv", 
#   target_glom = "Genus", 
#   significance_level = 0.05
# )
# 
# print(results_DESeq2$plot)
# head(results_DESeq2$results)  # View significant taxa
# results_DESeq2$ps_significant
# 
