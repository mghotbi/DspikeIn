#' @title Perform and Visualize Differential Abundance Analysis with edgeR or DESeq2
#' @description Performs differential abundance analysis using edgeR or DESeq2, processes data, and creates a volcano plot.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param method A string: `"edgeR"` or `"DESeq2"`.
#' @param group_var A string: grouping variable in the sample metadata.
#' @param contrast A character vector: levels to compare (e.g., `c("Control", "Treated")`).
#' @param pseudocount A numeric: pseudocount for zero handling (default = `1`).
#' @param significance_level A numeric: FDR threshold (default = `0.05`).
#' @param output_csv_path A string: output path for results CSV (optional).
#' @param target_glom A string: taxonomic rank to aggregate taxa (default = `"Genus"`).
#' @param palette A vector: colors for significant and non-significant points (default = `c("#FFEB3B", "#073B4C")`).
#'
#' @return A list with:
#'   \item{results}{Data frame of differential abundance results.}
#'   \item{obj_significant}{Filtered `phyloseq` or `TreeSummarizedExperiment` object.}
#'   \item{plot}{`ggplot2` volcano plot with factor levels.}
#'
#' @importFrom DspikeIn convert_tse_to_phyloseq convert_phyloseq_to_tse
#' @importFrom phyloseq otu_table sample_data tax_table prune_taxa tax_glom transform_sample_counts
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmFit glmLRT topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr mutate rename left_join filter select
#' @importFrom microbiome meta
#' @importFrom limma makeContrasts
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal scale_color_manual scale_shape_manual scale_size_continuous
#' @importFrom utils write.csv
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   results_edgeR <- perform_and_visualize_DA(
#'     obj = physeq_16SOTU,
#'     method = "edgeR",
#'     group_var = "Diet",
#'     contrast = c("Insectivore", "Carnivore"),
#'     output_csv_path = "DA_edgeR.csv",
#'     target_glom = "Genus",
#'     significance_level = 0.05
#'   )
#'
#'   # Print and visualize results
#'   print(results_edgeR$plot)
#'   head(results_edgeR$results)
#'   print(results_edgeR$obj_significant)
#'
#'   # Example with a TreeSummarizedExperiment (TSE) object
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   results_DESeq2 <- perform_and_visualize_DA(
#'     obj = tse_16SOTU,
#'     method = "DESeq2",
#'     group_var = "Diet",
#'     contrast = c("Insectivore", "Carnivore"),
#'     output_csv_path = "DA_DESeq2.csv",
#'     target_glom = "Genus",
#'     significance_level = 0.05
#'   )
#'
#'   # Print and visualize results
#'   print(results_DESeq2$plot)
#'   head(results_DESeq2$results)
#'   print(results_DESeq2$obj_significant)
#' }
#' }
#' @export
perform_and_visualize_DA <- function(obj, method, group_var, contrast,
                                     pseudocount = 1, significance_level = 0.05,
                                     output_csv_path = NULL,
                                     target_glom = "Genus",
                                     palette = c("#FFEB3B", "#073B4C")) {

  #  Detect if input is TSE and convert to phyloseq
  is_TSE <- inherits(obj, "TreeSummarizedExperiment")
  if (is_TSE) {
    obj <- DspikeIn::convert_tse_to_phyloseq(obj)
  }

  #  Extract and process Metadata
  metadata <- as.data.frame(microbiome::meta(obj))

  #  Ensure `group_var` exists*
  if (!group_var %in% colnames(metadata)) {
    stop(paste("Error: group_var", group_var, "not found in sample metadata!"))
  }

  # Convert "-" to "_" in metadata to comply with DESeq2
  metadata[[group_var]] <- factor(metadata[[group_var]])
  levels(metadata[[group_var]]) <- make.names(gsub("-", "_", levels(metadata[[group_var]])))

  # Apply fixed metadata back to `phyloseq` object
  phyloseq::sample_data(obj) <- phyloseq::sample_data(metadata)

  #  Remove zero/negative counts
  obj <- phyloseq::transform_sample_counts(obj, function(x) ifelse(x <= 0, pseudocount, x))

  #  Aggregate taxa at the specified rank
  obj <- phyloseq::tax_glom(obj, taxrank = target_glom)

  #  Extract Taxonomy Table
  taxonomy_table <- as.data.frame(phyloseq::tax_table(obj))
  taxonomy_table$OTU <- rownames(taxonomy_table)

  #  Convert contrast levels to valid names
  contrast_fixed <- make.names(gsub("-", "_", contrast))

  #  Check if contrast levels exist
  available_levels <- levels(metadata[[group_var]])
  if (!all(contrast_fixed %in% available_levels)) {
    stop(paste0(
      "Error: Specified contrast levels do not match group_var levels in metadata!\n",
      "Provided contrast: ", paste(contrast_fixed, collapse = ", "), "\n",
      "Available levels in metadata: ", paste(available_levels, collapse = ", "), "\n",
      "Ensure contrast matches the exact spelling in metadata."
    ))
  }

  ### ** edgeR Analysis**
  perform_edgeR <- function(obj, group_var, contrast_fixed, significance_level, taxonomy_table) {
    otu <- as(phyloseq::otu_table(obj), "matrix")

    metadata <- as.data.frame(microbiome::meta(obj))
    metadata[[group_var]] <- as.factor(metadata[[group_var]])

    dge <- edgeR::DGEList(counts = otu, group = metadata[[group_var]])
    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateDisp(dge)

    design <- model.matrix(~ 0 + metadata[[group_var]])
    colnames(design) <- levels(metadata[[group_var]])

    fit <- edgeR::glmFit(dge, design)
    contrast_matrix <- limma::makeContrasts(contrasts = paste0(contrast_fixed[2], "-", contrast_fixed[1]), levels = design)

    lrt <- edgeR::glmLRT(fit, contrast = contrast_matrix)
    res <- edgeR::topTags(lrt, n = Inf)$table %>%
      dplyr::rename(pvalue = PValue) %>%
      dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"),
                    Significance = ifelse(FDR < significance_level, "Significant", "Not Significant"),
                    group = ifelse(logFC > 0, contrast[2], contrast[1]))

    res$OTU <- rownames(res)
    res <- dplyr::left_join(res, taxonomy_table, by = "OTU")

    return(res)
  }

  ### ** DESeq2 Analysis**
  perform_DESeq2 <- function(obj, group_var, contrast_fixed, significance_level, taxonomy_table) {
    otu <- round(as(phyloseq::otu_table(obj), "matrix"))
    metadata <- as.data.frame(microbiome::meta(obj))

    #  Ensure metadata matches OTU table
    sample_ids <- sample_names(obj)
    metadata <- metadata[rownames(metadata) %in% sample_ids, , drop = FALSE]
    metadata <- metadata[match(sample_ids, rownames(metadata)), , drop = FALSE]

    #  Create DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, design = as.formula(paste("~", group_var)))
    dds <- DESeq2::DESeq(dds)

    #  Extract DE results
    res <- as.data.frame(DESeq2::results(dds, contrast = c(group_var, contrast_fixed[1], contrast_fixed[2])))

    res <- res %>%
      dplyr::rename(logFC = log2FoldChange, pvalue = pvalue) %>%
      dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"),
                    Significance = ifelse(FDR < significance_level, "Significant", "Not Significant"),
                    group = ifelse(logFC > 0, contrast[2], contrast[1]))

    res$OTU <- rownames(res)
    res <- dplyr::left_join(res, taxonomy_table, by = "OTU")

    return(res)
  }

  #  Perform DA Analysis
  results <- if (method == "edgeR") {
    perform_edgeR(obj, group_var, contrast_fixed, significance_level, taxonomy_table)
  } else {
    perform_DESeq2(obj, group_var, contrast_fixed, significance_level, taxonomy_table)
  }

  #  Remove NA & Apply Factor Levels
  results <- results %>%
    dplyr::filter(!is.na(Significance)) %>%
    dplyr::mutate(group = factor(group, levels = unique(group)))

  # Prune significant taxa
  obj_significant <- phyloseq::prune_taxa(results$OTU[results$Significance == "Significant"], obj)

  #  Volcano Plot
  p <- ggplot2::ggplot(results, ggplot2::aes(
    x = logFC,
    y = -log10(pvalue),
    color = Significance,
    shape = factor(group),
    size = -log10(FDR)
  )) +
    ggplot2::geom_point(alpha = 0.8) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = ifelse(FDR < significance_level & -log10(pvalue) > 5, OTU, "")),
      size = 4, fontface = "plain", max.overlaps = 15, segment.color = "grey40"
    ) +
    ggplot2::geom_hline(yintercept = -log10(significance_level),
                        linetype = "dashed", color = "#FF5733", linewidth = 1) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.text = ggplot2::element_text(size = 14, face = "plain"),
      axis.title = ggplot2::element_text(size = 14, face = "plain"),
      legend.position = "right",
      legend.box = "vertical",
      legend.key.size = ggplot2::unit(0.9, "cm"),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12, face = "plain"),
      legend.spacing.x = ggplot2::unit(-0.3, 'cm'),
      legend.spacing.y = ggplot2::unit(-0.1, 'cm')
        ) +
    ggplot2::labs(
      x = "Log2 Fold Change",
      y = "-log10 (P-value)",
      title = paste("Volcano Plot -", group_var),
      subtitle = paste("FDR threshold:", significance_level)
    )

  #  Convert back to TSE if needed
  if (is_TSE) {
    obj_significant <- DspikeIn::convert_phyloseq_to_tse(obj_significant)
  }

  return(list(results = results, obj_significant = obj_significant, plot = p))
}

# # Usage Example
# results_DESeq2 <- perform_and_visualize_DA(
#   obj = physeq_ITSOTU,
#   method = "DESeq2",
#   group_var = "Diet",
#   contrast = c("Insectivore", "Carnivore"),
#   output_csv_path = "DA_DESeq2.csv",
#   target_glom = "Genus",
#   significance_level = 0.05
# )

# print(results_DESeq2$plot)
# head(results_DESeq2$results)  # View significant taxa
# results_DESeq2$obj_significant

# perform_and_visualize_DA(
#   obj = tse_ITSOTU,
#   method = "DESeq2",
#   group_var = "Habitat",
#   contrast = c("Permanent Water", "Rainforest"),
#   significance_level = 0.05,
#   output_csv_path = "DA_DESeq2_results.csv",
#   target_glom = "Genus",
#   palette = DspikeIn::color_palette$MG  # Customizable colors
# )

# Run Differential Abundance Analysis
# results_edgeR <- perform_and_visualize_DA(
#   obj = tse_16SOTU,
#   method = "edgeR",
#   group_var = "Diet",
#   contrast = c("Insectivore", "Carnivore"),
#   output_csv_path = "DA_edgeR.csv",
#   target_glom = "Genus",
#   significance_level = 0.05
# )

# Display results
# print(results_edgeR$plot)
# head(results_edgeR$results)  # View significant taxa
# results_edgeR$obj_significant

# results_DESeq2 <- perform_and_visualize_DA(
#   obj = physeq_16SOTU,
#   method = "DESeq2",
#   group_var = "Diet",
#   contrast = c("Insectivore", "Carnivore"),
#   output_csv_path = "DA_DESeq2.csv",
#   target_glom = "Genus",
#   significance_level = 0.05
# )

# print(results_DESeq2$plot)
# head(results_DESeq2$results)  # View significant taxa
# results_DESeq2$obj_significant
