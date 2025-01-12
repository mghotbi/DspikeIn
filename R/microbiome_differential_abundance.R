#' Perform and Visualize Differential Abundance Analysis with edgeR or DESeq2
#'
#' This function performs differential abundance analysis using either edgeR or DESeq2, preprocesses the data,
#' and generates a volcano plot based on microbiome data.
#'
#' @param ps A `phyloseq` object containing microbiome data, including an OTU table and sample metadata.
#' @param method A character string specifying the method to use: \code{"edgeR"} or \code{"DESeq2"}.
#' @param group_var A string specifying the grouping variable in the sample metadata of the `phyloseq` object.
#' @param contrast A character vector specifying the levels to compare in the grouping variable (e.g., \code{c("Control", "Treated")}).
#' @param pseudocount A numeric value for pseudocount addition to handle zeros in count data. Default is \code{1}.
#' @param significance_level A numeric value specifying the significance threshold for adjusted p-values (FDR). Default is \code{0.05}.
#' @param output_csv_path A string specifying the path to save the results as a CSV file. Default is \code{"DA_results.csv"}.
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is \code{3}.
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa for analysis. Default is \code{"Genus"}.
#' @param palette A character vector of color hex codes for plotting. Default is \code{color_palette$extended_palette}.
#' 
#' @return A list with three components:
#' \itemize{
#'   \item \code{results}: A data frame containing differential abundance analysis results, including \code{logFC}, \code{pvalue}, \code{FDR}, \code{diff_abn}, and \code{treatment}.
#'   \item \code{ps_significant}: A `phyloseq` object containing only the significant OTUs/ASVs.
#'   \item \code{plot}: A ggplot2 object representing the volcano plot of the differential abundance results.
#' }
#'
#' @details
#' This function preprocesses microbiome count data, performs differential abundance analysis using 
#' edgeR or DESeq2, and visualizes the results as a volcano plot. Significant OTUs/ASVs are filtered 
#' based on the specified FDR threshold and included in the returned phyloseq object.
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
#' # Example with edgeR
#' results_edgeR <- perform_and_visualize_DA(
#'   ps = ps,
#'   method = "edgeR",
#'   group_var = "Treatment",
#'   contrast = c("Control", "Treated"),
#'   output_csv_path = "DA_edgeR.csv",
#'   target_glom = "Genus",
#'   significance_level = 0.05
#' )
#' print(results_edgeR$plot)
#' 
#' # Example with DESeq2
#' results_DESeq2 <- perform_and_visualize_DA(
#'   ps = ps,
#'   method = "DESeq2",
#'   group_var = "Treatment",
#'   contrast = c("Control", "Treated"),
#'   output_csv_path = "DA_DESeq2.csv",
#'   target_glom = "Genus",
#'   significance_level = 0.05
#' )
#' print(results_DESeq2$plot)
#' 
#' # Access significant taxa
#' significant_ps <- results_edgeR$ps_significant
#' head(results_edgeR$results)  # View significant taxa
#' }
#' 
#' @export 
perform_and_visualize_DA <- function(ps, method, group_var, contrast, pseudocount = 1, significance_level = 0.05, 
                                     output_csv_path = "DA_results.csv", point_size = 3, 
                                     target_glom = "Genus", palette = color_palette$extended_palette) {

remove_zero_negative_count_samples <- function(ps, pseudocount = 1) {
    # Extract the OTU table as a matrix/phyloseq
    otu <- as(phyloseq::otu_table(ps), "matrix")
    
    # Identify samples with zero or negative sums, or those containing NA values
    zero_negative_samples <- phyloseq::sample_sums(ps) <= 0 | apply(otu, 2, function(x) any(is.na(x)))
    
    # Remove problematic samples
    ps <- phyloseq::prune_samples(!zero_negative_samples, ps)
    
    # Add pseudocount to the OTU table
    otu <- as(phyloseq::otu_table(ps), "matrix") + pseudocount
    
    # Update the OTU table in the phyloseq object
    phyloseq::otu_table(ps) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
    
    # Return the cleaned phyloseq object
    return(ps)
  }
  
convert_categorical_to_factors <- function(ps) {
    metadata <- as(phyloseq::sample_data(ps), "data.frame")
    metadata[] <- lapply(metadata, function(x) if (is.character(x)) as.factor(x) else x)
    phyloseq::sample_data(ps) <- phyloseq::sample_data(metadata)
    return(ps)
  }
  
glom_taxa_at_rank <- function(physeq, rank = "Genus") {
    if (is.null(phyloseq::tax_table(physeq))) stop("A taxonomy table is required for tax_glom.")
    return(phyloseq::tax_glom(physeq, taxrank = rank))
  }
  
perform_DESeq2 <- function(ps, group_var, contrast, significance_level = 0.05) {
    otu <- as(phyloseq::otu_table(ps), "matrix")
    metadata <- as(phyloseq::sample_data(ps), "data.frame")
    metadata[[group_var]] <- as.factor(metadata[[group_var]])
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, design = as.formula(paste("~", group_var)))
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds, contrast = c(group_var, contrast[1], contrast[2]))
    res <- as.data.frame(res) %>%
      dplyr::rename(logFC = log2FoldChange, pvalue = pvalue) %>%
      dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"), diff_abn = FDR < significance_level) %>%
      dplyr::filter(!is.na(logFC) & !is.na(pvalue) & !is.na(padj))
    res$OTU <- rownames(res)
    return(res)
  }
  
perform_edgeR <- function(ps, group_var, contrast, significance_level = 0.05) {
    otu <- as(phyloseq::otu_table(ps), "matrix")
    metadata <- as(phyloseq::sample_data(ps), "data.frame")
    metadata[[group_var]] <- as.factor(metadata[[group_var]])
    dge <- edgeR::DGEList(counts = otu, group = metadata[[group_var]])
    design <- stats::model.matrix(~0 + metadata[[group_var]])
    colnames(design) <- make.names(levels(metadata[[group_var]]))
    contrast_matrix <- limma::makeContrasts(contrasts = paste(make.names(contrast[1]), "-", make.names(contrast[2])), 
                                            levels = colnames(design))
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmFit(dge, design)
    lrt <- edgeR::glmLRT(fit, contrast = contrast_matrix)
    res <- edgeR::topTags(lrt, n = Inf)$table %>%
      dplyr::rename(pvalue = PValue) %>%
      dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"), diff_abn = FDR < significance_level) %>%
      dplyr::filter(!is.na(logFC) & !is.na(pvalue) & !is.na(FDR))
    res$OTU <- rownames(res)
    return(res)
  }
  
  # Preprocess the phyloseq object
  ps <- remove_zero_negative_count_samples(ps, pseudocount = pseudocount)
  ps <- convert_categorical_to_factors(ps)
  ps <- glom_taxa_at_rank(ps, rank = target_glom)
  
  # Perform differential abundance analysis
  results <- if (method == "edgeR") {
    perform_edgeR(ps, group_var, contrast, significance_level)
  } else {
    perform_DESeq2(ps, group_var, contrast, significance_level)
  }
  
  # Save results to CSV
  utils::write.csv(results, output_csv_path, row.names = FALSE)
  
  # Filter significant OTUs
  significant_otus <- results %>% dplyr::filter(FDR < significance_level) %>% dplyr::pull(OTU)
  ps_significant <- phyloseq::prune_taxa(significant_otus, ps)
  
  # Merge results with taxonomy metadata
  tax_metadata <- as.data.frame(phyloseq::tax_table(ps_significant))
  tax_metadata$OTU <- rownames(tax_metadata)
  merged_results <- dplyr::left_join(results, tax_metadata, by = "OTU")
  
  # Add treatment column based on contrast
  merged_results <- merged_results %>%
    dplyr::mutate(treatment = ifelse(logFC > 0, contrast[2], contrast[1]))
  
  # Ensure pvalue is numeric and valid for plotting
  merged_results <- merged_results %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::mutate(pvalue = ifelse(is.numeric(pvalue) & pvalue <= 0, 1e-10, as.numeric(pvalue)))
  
  # Generate volcano plot
  p <- ggplot2::ggplot(merged_results, ggplot2::aes(x = logFC, y = -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(color = diff_abn, shape = treatment), size = point_size) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = paste("Volcano Plot -", group_var)) +
    ggplot2::theme(legend.position = "bottom")
  
  # Return results
  return(list(results = merged_results, ps_significant = ps_significant, plot = p))
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
