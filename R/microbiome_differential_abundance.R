#' Perform and Visualize Differential Abundance Analysis with edgeR
#'
#' This function performs edgeR analysis and generates a volcano plot based on microbiome data.
#'
#' @param ps A phyloseq object containing microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param contrast A vector specifying the levels to compare (e.g., c("Acris", "Anaxyrus")).
#' @param pseudocount A numeric value for pseudocount addition. Default is 1.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @param output_csv_path Path to save the edgeR results as a CSV file. Default is "DA_edgeR.csv".
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 3.
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa. Default is "Genus".
#' @param palette A character vector of color hex codes for plotting. Default is `MG()`.
#' @return A list containing the edgeR results, a phyloseq object with significant OTUs/ASVs, and a volcano plot.
#'
#' @importFrom phyloseq otu_table sample_data tax_table prune_taxa tax_glom
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom limma makeContrasts
#' @importFrom dplyr arrange mutate filter pull left_join rename
#' @importFrom ggplot2 ggplot aes geom_point scale_y_continuous theme_minimal labs scale_color_manual scale_shape_manual theme
#' @importFrom scales label_scientific
#' @importFrom utils write.csv 
#'
#' @examples
#' \dontrun{
#' results_edgeR <- perform_and_visualize_DA_edgeR(
#'   ps = ps, 
#'   group_var = "Host.genus", 
#'   contrast = c("Acris", "Anaxyrus"), 
#'   output_csv_path = "DA_edgeR.csv", 
#'   target_glom = "Genus", 
#'   significance_level = 0.05
#' )
#' 
#' # Visualize Volcano Plot
#' print(results_edgeR$plot)
#' }
#' 
#' @export
perform_and_visualize_DA_edgeR <- function(ps, group_var, contrast, pseudocount = 1, significance_level = 0.05, 
                                           output_csv_path = "DA_edgeR.csv", point_size = 3, 
                                           target_glom = "Genus", palette = MG()) {
  
  # Step 1: Remove samples with zero, negative counts, or NA values and add pseudocount
  ps <- remove_zero_negative_count_samples(ps, pseudocount = pseudocount)
  
  # Step 2: Convert categorical columns in sample data to factors
  ps <- convert_categorical_to_factors(ps)
  
  # Step 3: Glom taxa at the desired rank (e.g., Genus)
  ps <- glom_taxa_at_rank(ps, rank = target_glom)
  
  # Step 4: Ensure that the contrast levels are valid within the group_var
  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  if (!all(contrast %in% levels(metadata[[group_var]]))) {
    stop("One or both of the levels in the contrast do not exist in the group_var factor levels.")
  }
  
  # Step 5: Perform edgeR analysis
  results <- perform_edgeR(ps, group_var, contrast, significance_level = significance_level)
  
  # Step 6: Ensure that results contain logFC, pvalue, and FDR columns
  required_columns <- c("logFC", "PValue", "FDR")
  if (!all(required_columns %in% colnames(results))) {
    stop("Merged results are missing required columns (logFC, pvalue, FDR).")
  }
  
  # Step 7: Save results to CSV file
  utils::write.csv(results, output_csv_path, row.names = FALSE)
  cat("edgeR results saved to:", output_csv_path, "\n")
  
  # Step 8: Filter significant OTUs and create a phyloseq object with significant taxa
  significant_otus <- results %>% dplyr::filter(FDR < significance_level) %>% dplyr::pull(OTU)
  ps_significant <- phyloseq::prune_taxa(significant_otus, ps)
  
  # Step 9: Merge differential abundance results with phyloseq metadata
  tax_metadata <- as.data.frame(phyloseq::tax_table(ps_significant))
  tax_metadata$OTU <- rownames(tax_metadata)
  merged_results <- dplyr::left_join(results, tax_metadata, by = "OTU")
  
  # Step 10: Ensure columns are available after merging
  if (!all(required_columns %in% colnames(merged_results))) {
    stop("Merged results after join are missing required columns (logFC, pvalue, FDR).")
  }
  
  # Step 11: Add a column for the treatment level (shape by treatment)
  merged_results$treatment <- ifelse(merged_results$logFC > 0, contrast[1], contrast[2])
  
  # Step 12: Generate volcano plot
  p <- ggplot2::ggplot(merged_results, ggplot2::aes(x = logFC, y = -log10(PValue))) +
    ggplot2::geom_point(ggplot2::aes(color = diff_abn, shape = treatment), size = point_size) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +  # Shapes for treatment levels
    ggplot2::scale_color_manual(values = palette) +  # Use the MG color palette
    ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = paste("Volcano Plot -", group_var)) +
    ggplot2::theme(legend.position = "bottom")
  
  return(list(results = merged_results, ps_significant = ps_significant, plot = p))
}

#' Perform and Visualize Differential Abundance Analysis with DESeq2
#'
#' This function performs DESeq2 analysis and generates a volcano plot based on microbiome data.
#'
#' @param ps A phyloseq object containing microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param contrast A vector specifying the levels to compare (e.g., c("Acris", "Anaxyrus")).
#' @param pseudocount A numeric value for pseudocount addition. Default is 1.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @param output_csv_path Path to save the DESeq2 results as a CSV file. Default is "DA_deseq2.csv".
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 3.
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa. Default is "Genus".
#' @param palette A character vector of color hex codes for plotting. Default is `MG()`.
#' @return A list containing the DESeq2 results, a phyloseq object with significant OTUs/ASVs, and a volcano plot.
#'
#' @importFrom phyloseq otu_table sample_data tax_table prune_taxa tax_glom
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr filter pull left_join rename
#' @importFrom ggplot2 ggplot aes geom_point scale_y_continuous theme_minimal labs scale_color_manual scale_shape_manual
#' @importFrom scales label_scientific
#' @importFrom utils write.csv 
#' @examples
#' \dontrun{
#' results_DESeq2 <- perform_and_visualize_DA_DESeq2(
#' ps = ps, 
#' group_var = "Host.genus", 
#' contrast = c("Acris", "Anaxyrus"), 
#' output_csv_path = "DA_deseq2.csv", 
#' target_glom = "Genus", 
#' significance_level = 0.05
#' )
#' # Visualize Volcano Plot
#' print(results_DESeq2$plot)
#' results_deseq2$results
#' results_deseq2$ps_significant
#' }
#' @export
perform_and_visualize_DA_DESeq2 <- function(ps, group_var, contrast, pseudocount = 1, significance_level = 0.05, 
                                            output_csv_path = "DA_deseq2.csv", point_size = 3, 
                                            target_glom = "Genus", palette = MG()) {
  
  # Step 1: Remove samples with zero, negative counts, or NA values and add pseudocount
  ps <- remove_zero_negative_count_samples(ps, pseudocount = pseudocount)
  
  # Step 2: Convert categorical columns in sample data to factors
  ps <- convert_categorical_to_factors(ps)
  
  # Step 3: Glom taxa at the desired rank (e.g., Genus)
  ps <- glom_taxa_at_rank(ps, rank = target_glom)
  
  # Step 4: Perform DESeq2 analysis
  results <- perform_DESeq2(ps, group_var, contrast, significance_level = significance_level)
  
  # Step 5: Save results to CSV file
  utils::write.csv(results, output_csv_path, row.names = FALSE)
  cat("DESeq2 results saved to:", output_csv_path, "\n")
  
  # Step 6: Filter significant OTUs and create a phyloseq object with significant taxa
  significant_otus <- results %>% dplyr::filter(FDR < significance_level) %>% dplyr::pull(OTU)
  ps_significant <- phyloseq::prune_taxa(significant_otus, ps)
  
  # Step 7: Merge differential abundance results with phyloseq metadata
  tax_metadata <- as.data.frame(phyloseq::tax_table(ps_significant))
  tax_metadata$OTU <- rownames(tax_metadata)
  merged_results <- dplyr::left_join(results, tax_metadata, by = "OTU")
  
  # Step 8: Ensure columns are available after merging
  required_columns <- c("logFC", "pvalue", "FDR")
  if (!all(required_columns %in% colnames(merged_results))) {
    stop("Merged results are missing required columns (logFC, pvalue, FDR).")
  }
  
  # Step 9: Add a column for the treatment level (shape by treatment)
  merged_results$treatment <- ifelse(merged_results$logFC > 0, contrast[1], contrast[2])
  
  # Step 10: Generate volcano plot
  p <- ggplot2::ggplot(merged_results, ggplot2::aes(x = logFC, y = -log10(pvalue))) +
    ggplot2::geom_point(ggplot2::aes(color = FDR < significance_level, shape = treatment), size = point_size) +
    ggplot2::scale_shape_manual(values = c(16, 17)) +  # Shapes for treatment levels
    ggplot2::scale_color_manual(values = palette) +  # Color based on significance
    ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Log2 Fold Change", y = "-log10(P-value)", title = paste("Volcano Plot -", group_var)) +
    ggplot2::theme(legend.position = "bottom")
  
  return(list(results = merged_results, ps_significant = ps_significant, plot = p))
}


#' Perform Pairwise Comparisons with DESeq2
#'
#' This function performs pairwise comparisons for differential abundance analysis using DESeq2.
#' It filters out NA values from the results before returning the data frame.
#'
#' @param ps A phyloseq object containing microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param contrast A vector specifying the levels to compare (e.g., c("level1", "level2")).
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @return A data frame containing the differential abundance results.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom dplyr arrange mutate pull rename filter
#' @export
perform_DESeq2 <- function(ps, group_var, contrast, significance_level = 0.05) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  
  # Ensure grouping variable is a factor
  metadata[[group_var]] <- as.factor(metadata[[group_var]])
  
  # Validate that the levels in the contrast exist in the factor levels
  if (!all(contrast %in% levels(metadata[[group_var]]))) {
    stop("One or both of the levels in the contrast do not exist in the group_var factor levels.")
  }
  
  # Create DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, design = as.formula(paste("~", group_var)))
  
  # Run DESeq2
  dds <- DESeq2::DESeq(dds)
  
  # Perform custom contrast comparison
  res <- DESeq2::results(dds, contrast = c(group_var, contrast[1], contrast[2]))
  res <- as.data.frame(res)
  
  # Rename log2FoldChange to logFC for consistency with edgeR
  res <- res %>%
    dplyr::rename(logFC = log2FoldChange)
  
  # Add OTU column and remove rows with NA in critical columns
  res$OTU <- rownames(res)
  res <- res %>%
    dplyr::filter(!is.na(logFC) & !is.na(pvalue) & !is.na(padj)) %>%  # Filter out rows with NA values
    dplyr::mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    dplyr::mutate(diff_abn = FDR < significance_level)
  
  return(res)
}

#' Perform Pairwise Comparisons with edgeR
#'
#' This function performs pairwise comparisons for differential abundance analysis using edgeR.
#' It filters out NA values from the results before returning the data frame.
#'
#' @param ps A phyloseq object containing microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param contrast A vector specifying the levels to compare (e.g., c("level1", "level2")).
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @return A data frame containing the differential abundance results.
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
#' @importFrom limma makeContrasts
#' @importFrom dplyr arrange mutate pull rename filter
#' @export
perform_edgeR <- function(ps, group_var, contrast, significance_level = 0.05) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  metadata <- as(phyloseq::sample_data(ps), "data.frame")
  
  # Ensure grouping variable is a factor
  metadata[[group_var]] <- as.factor(metadata[[group_var]])
  
  # Check that contrast levels exist in the group_var
  if (!all(contrast %in% levels(metadata[[group_var]]))) {
    stop("One or both of the levels in the contrast do not exist in the group_var factor levels.")
  }
  
  # Create the DGEList object
  dge <- edgeR::DGEList(counts = otu, group = metadata[[group_var]])
  
  # Create the design matrix
  design <- stats::model.matrix(~ 0 + metadata[[group_var]])
  colnames(design) <- make.names(levels(metadata[[group_var]]))  # Ensure syntactically valid names
  
  # Print the column names to debug
  cat("Column names in the design matrix:\n", colnames(design), "\n")
  
  # Ensure contrast levels use syntactically valid names
  contrast1 <- make.names(contrast[1])
  contrast2 <- make.names(contrast[2])
  
  # Check if the contrast levels exist in the design matrix
  if (!contrast1 %in% colnames(design) || !contrast2 %in% colnames(design)) {
    stop(paste("One or both contrasts", contrast1, "or", contrast2, "do not exist in the design matrix."))
  }
  
  # Create the contrast expression
  contrast_expr <- paste0(contrast1, "-", contrast2)
  
  # Generate the contrast matrix
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_expr, levels = colnames(design))
  
  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit the model using glmFit
  fit <- edgeR::glmFit(dge, design)
  
  # Perform the likelihood ratio test using the contrast matrix
  lrt <- edgeR::glmLRT(fit, contrast = contrast_matrix)
  
  # Extract the top results
  res <- edgeR::topTags(lrt, n = Inf)$table
  res$OTU <- rownames(res)
  res$comparison <- paste0(contrast[1], "_vs_", contrast[2])
  
  res <- res %>%
    dplyr::arrange(PValue) %>%
    dplyr::mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    dplyr::mutate(diff_abn = FDR < significance_level)
  
  # Remove rows with NA values in logFC, PValue, or FDR
  res <- res %>%
    dplyr::filter(!is.na(logFC) & !is.na(PValue) & !is.na(FDR))
  
  return(res)
}

# Supporting functions -----------------------------------------------------
#' 
#' MG Color Palette
#' This function returns a character vector of color palettes used in the package.
#' @return A character vector of color hex codes.
#' @export
MG <- function() {
  c(
    "#FFFF33", "#FF7F00", "#E41A1C", "firebrick4", "#2e4057", "#984EA3", "#377EB8", "olivedrab3",
    "#4DAF4A", "#336633", "grey80", "#BB650B", "gold", "#559999", "#7570b3", "#E78AC3", "#A6D854",
    "#66a61e", "#e6ab02", "#a6761d", "#663300", "#66C2A5", "#0e669b", "#00798c", "dodgerblue4",
    "steelblue2", "#00AFBB", "#E7B800", "#FC4E07", "lightskyblue4", "green", "red", "#FFF000",
    "#0099CC", "#FF9933", "#CC9900", "chartreuse1", "#FF3399", "#00FFFF", "#0000CC", "#A37F6F",
    "#9183E6", "#00AD9A", "#990033", "#909800", "#00FF00", "#17b5b4", "#AED1D6", "#b1010c",
    "firebrick2", "blue", "navy", "yellow", "brown", "black", "purple", "darkred", "darkgreen",
    "#82cfd0", "#b2e0e4", "honeydew3", "#8d96a3", "lavender", "#CC6686", "lavenderblush2",
    "mistyrose3", "#e1deda", "darkgoldenrod", "burlywood", "papayawhip", "wheat4", "cornsilk3",
    "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream",
    "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue"  
  )
}

#' Remove Samples with Zero, Negative Counts, or NA Values and Add Pseudocount
#'
#' This function removes samples with zero, negative counts, or NA values and adds pseudocounts.
#' 
#' @param ps A phyloseq object containing the OTU table and sample data.
#' @param pseudocount A numeric value to add to avoid zero counts. Default is 1.
#' @return A phyloseq object with updated OTU table after pseudocount addition.
#' @importFrom phyloseq otu_table prune_samples prune_taxa sample_sums
#' @export
remove_zero_negative_count_samples <- function(ps, pseudocount = 1) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  
  # Remove samples with zero or negative counts
  zero_negative_count_samples <- phyloseq::sample_sums(ps) <= 0
  na_count_samples <- apply(otu, 2, function(x) any(is.na(x)))
  samples_to_remove <- zero_negative_count_samples | na_count_samples
  if (any(samples_to_remove)) {
    cat("Removing", sum(samples_to_remove), "samples with zero, negative counts, or NA values.\n")
    ps <- phyloseq::prune_samples(!samples_to_remove, ps)
    otu <- as(phyloseq::otu_table(ps), "matrix")
  }
  
  # Remove rows (features) with zero counts across all samples
  zero_rows <- rowSums(otu) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples.\n")
    otu <- otu[!zero_rows, ]
    ps <- phyloseq::prune_taxa(!zero_rows, ps)
  }
  
  # Add pseudocount and round to integer
  otu <- otu + pseudocount
  otu <- round(otu)
  phyloseq::otu_table(ps) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  return(ps)
}

#' Convert Categorical Columns to Factors in Sample Data
#'
#' This function converts categorical columns in sample data to factors.
#'
#' @param ps A phyloseq object.
#' @return A phyloseq object with updated sample data.
#' @importFrom phyloseq sample_data
#' @export
convert_categorical_to_factors <- function(ps) {
  sample_data_df <- as(phyloseq::sample_data(ps), "data.frame")
  for (col in colnames(sample_data_df)) {
    if (is.character(sample_data_df[[col]]) || is.factor(sample_data_df[[col]])) {
      sample_data_df[[col]] <- as.factor(sample_data_df[[col]])
    }
  }
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sample_data_df)
  return(ps)
}

#' Glom Taxa at Specific Level
#'
#' This function gloms taxa at a specific taxonomic level in a phyloseq object.
#'
#' @param physeq A phyloseq object containing the OTU and taxonomy tables.
#' @param rank A character string specifying the taxonomic rank to glom. Default is "Genus".
#' @return A phyloseq object with taxa glommed at the specified rank.
#' @importFrom phyloseq tax_glom tax_table ntaxa
#' @export
glom_taxa_at_rank <- function(physeq, rank = "Genus") {
  if (is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
    stop("The tax_glom() function requires that physeq contains a taxonomy table.")
  }
  
  physeq <- phyloseq::tax_glom(physeq, taxrank = rank)
  if (phyloseq::ntaxa(physeq) == 0) stop("No taxa remain after glomming.")
  
  return(physeq)
}

# Usage Examples,
# Example DESeq2,
# results_deseq2 <- perform_and_visualize_DA_DESeq2(
#  ps = ps, 
#  group_var = "Host.genus", 
#  contrast = c("Acris", "Anaxyrus"), 
#  output_csv_path = "DA_deseq2.csv", 
# target_glom = "Genus", 
#  significance_level = 0.05
#)
# print(results_deseq2$plot)
# results_deseq2$results
# results_deseq2$ps_significant


# Example edgeR,

# group_var <- "Host.genus"     # sample variable
# contrast <- c("Acris", "Anaxyrus")  # Specify the levels to compare

# with edgeR and visualize results
# results_edgeR <- perform_and_visualize_DA_edgeR(
#  ps = ps,                           #  phyloseq object
#  group_var = group_var,              # Grouping variable
#  contrast = contrast,                # Contrast levels
#  output_csv_path = "DA_edgeR.csv",   # Path to save the edgeR results
# target_glom = "Genus",              # Taxonomic rank to aggregate
#  significance_level = 0.05           # Significance threshold (default is 0.05)
#)

#print(results_edgeR$plot)
#head(results_edgeR$results)           # View the differential abundance results
#results_edgeR$ps_significant          # View the significant taxa (phyloseq obj)
