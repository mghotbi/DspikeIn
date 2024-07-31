<<<<<<< HEAD
# Load required packages
#' Install and Load Required Packages
#'
=======
>>>>>>> 9c89fa128aa50dce648041926482a4901ad3019b
#' This function installs and loads required packages.
#'
#' @param packages A vector of package names to be installed and loaded.
#' @examples
#' required_packages <- c("phyloseq", "DESeq2", "edgeR", "sva", "EDASeq", "RUVSeq", "BiocManager", "BiocGenerics", "ANCOMBC", "ggplot2", "dplyr", "DT")
#' install_and_load(required_packages)
#' @export
install_and_load <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("phyloseq", "DESeq2", "edgeR", "sva", "EDASeq", "RUVSeq", "BiocManager", "BiocGenerics", "ANCOMBC")) {
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

required_packages <- c("phyloseq", "DESeq2", "edgeR", "sva", "EDASeq", "RUVSeq", "BiocManager", "BiocGenerics", "ANCOMBC", "ggplot2", "dplyr", "DT")
install_and_load(required_packages)

#' MG Color Palette
#'
#' This function returns a character vector of color palettes used in the package.
#'
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
    "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream", "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue"
  )
}

#' Remove Samples with Zero, Negative Counts, or NA Values and Add Pseudocount
#'
#' @param ps A phyloseq object.
#' @param pseudocount A numeric value to add to avoid zero counts.
#' @return A phyloseq object with filtered and adjusted OTU table.
#' @examples
#' ps <- remove_zero_negative_count_samples(ps)
#' @export
remove_zero_negative_count_samples <- function(ps, pseudocount = 1e-6) {
  otu <- as(otu_table(ps), "matrix")
  
  zero_negative_count_samples <- sample_sums(ps) <= 0
  na_count_samples <- apply(otu, 2, function(x) any(is.na(x)))
  samples_to_remove <- zero_negative_count_samples | na_count_samples
  if (any(samples_to_remove)) {
    cat("Removing", sum(samples_to_remove), "samples with zero, negative counts, or NA values.\n")
    ps <- prune_samples(!samples_to_remove, ps)
    otu <- as(otu_table(ps), "matrix")
  }
  
  zero_rows <- rowSums(otu) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples.\n")
    otu <- otu[!zero_rows, ]
    ps <- prune_taxa(!zero_rows, ps)
  }
  
  otu <- otu + pseudocount
  otu <- round(otu)
  
  otu_table(ps) <- otu_table(otu, taxa_are_rows = TRUE)
  
  return(ps)
}

#' Convert Categorical Columns to Factors in Sample Data
#'
#' This function converts all character and factor columns in the sample data of a phyloseq object to factors.
#'
#' @param ps A phyloseq object.
#' @return A phyloseq object with updated sample data.
#' @examples
#' ps <- convert_categorical_to_factors(ps)
#' @export
convert_categorical_to_factors <- function(ps) {
  sample_data_df <- as(sample_data(ps), "data.frame")
  for (col in colnames(sample_data_df)) {
    if (is.character(sample_data_df[[col]]) || is.factor(sample_data_df[[col]])) {
      sample_data_df[[col]] <- as.factor(sample_data_df[[col]])
    }
  }
  sample_data(ps) <- sample_data(sample_data_df)
  return(ps)
}

#' Relativized Filtered Taxa
#'
#' This function filters taxa from a phyloseq object based on custom thresholds for percentage of samples, mean abundance, count, and relative abundance.
#'
#' @param physeq A phyloseq object containing the microbial data.
#' @param threshold_percentage A numeric value specifying the minimum percentage of samples in which a taxon must be present to be retained. Default is 0.5.
#' @param threshold_mean_abundance A numeric value specifying the minimum mean abundance of a taxon to be retained. Default is 0.001.
#' @param threshold_count A numeric value specifying the minimum count of a taxon in a sample to be considered present. Default is 10.
#' @param threshold_relative_abundance A numeric value specifying the minimum relative abundance of a taxon to be retained. Default is NULL.
#' @return A phyloseq object containing only the taxa that meet the specified thresholds.
#' @examples
#' FT <- relativized_filtered_taxa(spiked_16S, threshold_percentage = 0.6, threshold_mean_abundance = 0.0005, threshold_count = 5, threshold_relative_abundance = 0.01)
#' @export
relativized_filtered_taxa <- function(physeq, threshold_percentage = 0.5, threshold_mean_abundance = 0.001, threshold_count = 10, threshold_relative_abundance = NULL) {
  
  nsamples <- nsamples(physeq)
  sample_sum <- sample_sums(physeq)
  
  filter_function <- function(x) {
    (sum(x > threshold_count) > nsamples * threshold_percentage) | 
      ((sum(x > threshold_count) > (nsamples * 0.1)) & (mean(x / sample_sum) > threshold_mean_abundance) & (max(x / sample_sum) > threshold_relative_abundance))
  }
  
  two_way_filtered <- filter_taxa(physeq, filter_function, prune = TRUE)
  return(two_way_filtered)
}

#' Glom Taxa at Specific Level
#'
#' This function gloms taxa at a specified taxonomic rank.
#'
#' @param physeq A phyloseq object.
#' @param rank A character string specifying the taxonomic rank to glom.
#' @return A phyloseq object with taxa glommed at the specified rank.
#' @examples
#' ps_glommed <- glom_taxa_at_rank(ps, "Genus")
#' @export
glom_taxa_at_rank <- function(physeq, rank) {
  physeq <- tax_glom(physeq, taxrank = rank)
  return(physeq)
}
#' Perform Differential Abundance Analysis with edgeR
#'
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
perform_edgeR <- function(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level) {
  otu <- as(otu_table(ps), "matrix")
  metadata <- as(sample_data(ps), "data.frame")
  
  if (!group_var %in% colnames(metadata)) {
    stop("The specified group variable does not exist in the sample data.")
  }
  
  ps_filtered <- tryCatch({
    relativized_filtered_taxa(ps, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance)
  }, error = function(e) {
    warning("Filtering resulted in no taxa remaining: ", e$message, ". Skipping filtering step.")
    return(ps)
  })
  
  if (nsamples(ps_filtered) == 0 || ntaxa(ps_filtered) == 0) {
    stop("No taxa remain after filtering. Adjust the filtering thresholds.")
  }
  
  # Glom taxa at the specified rank
  ps_filtered_glom <- glom_taxa_at_rank(ps_filtered, "Genus")
  
  otu_normalized <- as(otu_table(ps_filtered_glom), "matrix")
  group <- metadata[[group_var]]
  
  dge <- DGEList(counts = otu_normalized, group = group)
  design <- model.matrix(~group)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
  
  top_tags <- topTags(lrt, n = Inf)
  results <- top_tags$table
  results$OTU <- rownames(results)
  results <- results %>%
    arrange(PValue) %>%
    mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    mutate(diff_abn = FDR < significance_level)
  
  return(results)
}

#' Perform Differential Abundance Analysis with DESeq2
#'
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
perform_DESeq2 <- function(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level) {
  otu <- as(otu_table(ps), "matrix")
  metadata <- as(sample_data(ps), "data.frame")
  
  if (!group_var %in% colnames(metadata)) {
    stop("The specified group variable does not exist in the sample data.")
  }
  
  ps_filtered <- tryCatch({
    relativized_filtered_taxa(ps, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance)
  }, error = function(e) {
    warning("Filtering resulted in no taxa remaining: ", e$message, ". Skipping filtering step.")
    return(ps)
  })
  
  if (nsamples(ps_filtered) == 0 || ntaxa(ps_filtered) == 0) {
    stop("No taxa remain after filtering. Adjust the filtering thresholds.")
  }
  
  # Glom taxa at the specified rank
  ps_filtered_glom <- glom_taxa_at_rank(ps_filtered, "Genus")
  
  otu <- as(otu_table(ps_filtered_glom), "matrix")
  otu <- round(otu)
  otu_table(ps_filtered_glom) <- otu_table(otu, taxa_are_rows = TRUE)
  
  dds <- DESeqDataSetFromMatrix(countData = otu, 
                                colData = metadata, 
                                design = as.formula(paste("~", group_var)))
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj), ]
  res <- as.data.frame(res)
  res$OTU <- rownames(res)
  res <- res %>%
    mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue) %>%
    mutate(logFC = log2FoldChange, PValue = pvalue) %>%
    mutate(diff_abn = FDR < significance_level)
  return(res)
}

#' Perform Differential Abundance Analysis with ANCOMBC
#'
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
perform_ANCOMBC <- function(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level) {
  ps_filtered <- tryCatch({
    relativized_filtered_taxa(ps, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance)
  }, error = function(e) {
    warning("Filtering resulted in no taxa remaining: ", e$message, ". Skipping filtering step.")
    return(ps)
  })
  
  # Glom taxa at the specified rank
  ps_filtered_glom <- glom_taxa_at_rank(ps, "Genus")
  
  result <- ANCOMBC::ancombc(phyloseq = ps_filtered_glom, formula = as.formula(paste("~", group_var)), 
                             p_adj_method = "holm", lib_cut = 1000, group = group_var, 
                             struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, max_iter = 100, 
                             conserve = TRUE, alpha = 0.05, global = TRUE)
  significant_otus <- result$res$diff_abn %>% filter(diff_abn == TRUE) %>% pull(taxon)
  results <- data.frame(OTU = significant_otus, diff_abn = TRUE, stringsAsFactors = FALSE)
  return(results)
}

#' Perform Differential Abundance Analysis
#'
#' This function normalizes and filters the data, performs differential abundance analysis, extracts significant OTUs,
#' rebuilds the phyloseq object, and merges results with metadata.
#'
#' @param ps A phyloseq object containing the microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param method A string specifying the differential abundance method ("edgeR", "DESeq2", "ANCOMBC"). Default is "edgeR".
#' @param threshold_percentage A numeric value for filtering threshold percentage. Default is 0.0001.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance. Default is 0.00001.
#' @param threshold_count A numeric value for filtering count. Default is 1.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance. Default is 0.00001.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @return A list containing the final results and the phyloseq object with significant OTUs.
#' @examples
#' result <- perform_differential_abundance(ps, group_var = "Animal.type", method = "edgeR", threshold_percentage = 0.0001, threshold_mean_abundance = 0.00001, threshold_count = 1, threshold_relative_abundance = 0.00001, significance_level = 0.05)
#' @export
perform_differential_abundance <- function(ps, group_var, method = "edgeR",
                                           threshold_percentage = 0.0001, threshold_mean_abundance = 0.00001, 
                                           threshold_count = 1, threshold_relative_abundance = 0.00001,
                                           significance_level = 0.05) {
  if (method == "edgeR") {
    results <- perform_edgeR(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level)
  } else if (method == "DESeq2") {
    results <- perform_DESeq2(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level)
  } else if (method == "ANCOMBC") {
    results <- perform_ANCOMBC(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level)
  } else {
    stop("Invalid differential abundance method")
  }
  
  significant_otus <- results %>%
    filter(FDR < significance_level) %>%
    pull(OTU)
  
  ps_significant <- prune_taxa(significant_otus, ps)
  
  melted_ps <- psmelt(ps_significant)
  final_results <- merge(results, melted_ps, by = "OTU")
  
  return(list(final_results = final_results, ps_significant = ps_significant))
}

#' Visualize Differential Abundance Results
#'
#' @param final_results The final results from the differential abundance analysis.
#' @param group_var The grouping variable used in the analysis.
#' @param point_size The size of points in the volcano plot. Default is 2.
#' @param palette A character vector specifying the color palette. Default is MG().
#' @return A ggplot2 object representing the volcano plot.
#' @export
visualize_differential_abundance <- function(final_results, group_var, point_size = 2, palette = MG()) {
  ggplot(final_results, aes(x = logFC, y = -log10(PValue), color = diff_abn)) +
    geom_point(size = point_size) +
    theme_minimal() +
    labs(x = "Log Fold Change", y = "-Log10 P-Value") +
    scale_color_manual(values = palette)
}

#' Generate a Taxa Barplot
#'
#' This function creates a bar plot of the relative or absolute abundances of taxa at a specified taxonomic rank.
#' The top taxa are selected, and the plot can be customized with various options.
#' The function uses the `tax_glom` function from the `phyloseq` package for taxonomic grouping.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param target_glom A character string specifying the taxonomic rank to plot (e.g., "Genus").
#'        This specifies the taxonomic level at which the data should be aggregated and visualized.
#' @param custom_tax_names A character vector specifying custom taxonomic names for the levels. Default is NULL.
#'        This allows the user to provide custom names for the taxonomic levels in the taxonomy table of the phyloseq object.
#' @param normalize A logical indicating whether to normalize the sample counts to relative abundances. Default is TRUE.
#' @param treatment_variable A character string specifying the treatment variable to use for the x-axis. Default is "Treatment".
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "relative".
#' @param x_angle A numeric value specifying the angle of the x-axis text labels. Default is 20.
#' @param fill_variable A character string specifying the variable to use for fill in stacking taxa. Default is target_glom.
#' @param facet_variable A character string specifying the variable to use for faceting. Default is "Phylum".
#' @param top_n_taxa A numeric value specifying the number of top taxa to include in the plot. Default is 20.
#' @return A list containing the ggplot2 bar plot object (`barplot`) and the pruned phyloseq object (`taxa_data`) with the top taxa.
#' @examples
#' bp <- taxa_barplot(physeqASV16, target_glom = "Genus", normalize = TRUE, treatment_variable = "animal.type", abundance_type = "relative")
#' print(bp$barplot)
#' @export
taxa_barplot <- function(physeq, target_glom = "Genus", custom_tax_names = NULL, normalize = TRUE, treatment_variable = "Treatment", abundance_type = "relative", x_angle = 20, fill_variable = target_glom, facet_variable = "Phylum", top_n_taxa = 20, palette = MG()) {
  
  glom <- tax_glom(physeq, taxrank = target_glom)
  glom_1 <- prune_taxa(taxa_sums(glom) > 0, glom)
  glom_C <- tax_glom(glom_1, taxrank = target_glom)
  
  top_taxa <- names(sort(taxa_sums(glom_C), decreasing = TRUE)[1:top_n_taxa])
  top_taxa_pruned <- prune_taxa(top_taxa, glom_C)
  top_v5 <- prune_taxa(taxa_sums(top_taxa_pruned) > 0, top_taxa_pruned)
  
  if (!is.null(custom_tax_names)) {
    colnames(tax_table(top_v5)) <- custom_tax_names
  } else {
    colnames(tax_table(top_v5)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", target_glom)
  }
  
  if (normalize && abundance_type == "relative") {
    top_v5 <- transform_sample_counts(top_v5, function(OTU) OTU / sum(OTU))
  }
  
  pm <- psmelt(top_v5)
  
  if (abundance_type == "relative") {
    p <- ggplot(pm, aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format())
  } else {
    p <- ggplot(pm, aes_string(x = treatment_variable, y = "Abundance", fill = fill_variable)) +
      geom_bar(stat = "identity")
  }
  
  p <- p + 
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 20, keyheight = unit(0.6, "lines"), keywidth = unit(0.6, "lines"))) +
    theme_minimal() +
    labs(x = "", y = if (abundance_type == "relative") "Relative Abundance" else "Absolute Abundance") +
    facet_grid(cols = vars(.data[[facet_variable]]), scales = "free") +
    theme(
      strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = x_angle, vjust = 0.5, hjust = 1)
    ) +
    scale_fill_manual(values = palette)
  
  return(list(barplot = p, taxa_data = top_v5))
}

#' Perform and Visualize Differential Abundance Analysis
#'
#' This function normalizes and filters the data, performs differential abundance analysis, extracts significant OTUs,
#' rebuilds the phyloseq object, merges results with metadata, and visualizes the results.
#'
#' @param ps A phyloseq object containing the microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param method A string specifying the differential abundance method ("edgeR", "DESeq2", "ANCOMBC"). Default is "edgeR".
#' @param threshold_percentage A numeric value for filtering threshold percentage. Default is 0.0001.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance. Default is 0.00001.
#' @param threshold_count A numeric value for filtering count. Default is 1.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance. Default is 0.00001.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 2.
#' @return A list containing the final results, the phyloseq object with significant OTUs, and ggplot objects for the visualizations.
#' @examples
#' ps <- remove_zero_negative_count_samples(ps)
#' results <- perform_and_visualize_differential_abundance(ps, "Animal.ecomode", method = "DESeq2", threshold_percentage = 0.001, threshold_mean_abundance = 0.001, threshold_count = 1, threshold_relative_abundance = 0.001, significance_level = 0.05, point_size = 2)
#' print(results$plot)
#' print(results$barplot)
#' @export
perform_and_visualize_differential_abundance <- function(ps, group_var, method = "edgeR",
                                                         threshold_percentage = 0.0001, threshold_mean_abundance = 0.00001, 
                                                         threshold_count = 1, threshold_relative_abundance = 0.00001,
                                                         significance_level = 0.05, point_size = 2) {
  
  results <- perform_differential_abundance(ps, group_var, method, threshold_percentage, 
                                            threshold_mean_abundance, threshold_count, threshold_relative_abundance, 
                                            significance_level)
  
  p <- visualize_differential_abundance(results$final_results, group_var, point_size = point_size, palette = MG())
  
  bp_rel <- taxa_barplot(results$ps_significant, target_glom = "Genus", treatment_variable = group_var, 
                         abundance_type = "relative", x_angle = 90, fill_variable = "Genus", facet_variable = "Phylum", top_n_taxa = 20, palette = MG())
  
  bp_abs <- taxa_barplot(results$ps_significant, target_glom = "Genus", treatment_variable = group_var, 
                         abundance_type = "absolute", x_angle = 90, fill_variable = "Genus", facet_variable = "Phylum", top_n_taxa = 20, palette = MG())
  
  return(list(plot = p, barplot_rel = bp_rel$barplot, barplot_abs = bp_abs$barplot, results = results))
}

# Example usage:
# ps is a phyloseq obj

# ps <- remove_zero_negative_count_samples(ps)
# ps <- convert_categorical_to_factors(ps)

# #edgeR example
# results <- perform_and_visualize_differential_abundance(ps, "Result", method = "edgeR", 
#                                                         threshold_percentage = 0.001, threshold_mean_abundance = 0.001, 
#                                                         threshold_count = 5, threshold_relative_abundance = 0.001, 
#                                                         significance_level = 0.05, point_size = 3)
# print(results$plot)
# print(results$barplot_rel)
# print(results$barplot_abs)
# 
# # DESeq2 example
# results <- perform_and_visualize_differential_abundance(ps, "Host.genus", method = "DESeq2", 
#                                                         threshold_percentage = 0.001, threshold_mean_abundance = 0.001, 
#                                                         threshold_count = 5, threshold_relative_abundance = 0.001, 
#                                                         significance_level = 0.05, point_size = 3)
# print(results$plot)
# print(results$barplot_rel)
# print(results$barplot_abs)
