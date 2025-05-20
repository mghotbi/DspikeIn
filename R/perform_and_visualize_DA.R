#' @title Perform and Visualize Differential Abundance Analysis with edgeR or DESeq2
#'
#' @description
#' Performs differential abundance analysis using edgeR or DESeq2, and visualizes results
#' with a volcano plot, a log-fold change bar plot, and a relative abundance bar plot
#' for significantly enriched taxa. Supports comparisons across one or more variables and
#' allows global FDR correction across multiple contrasts.
#'
#' @details
#' For edgeR, the standard error of log-fold change (lfcSE) is estimated using the formula:
#' \code{lfcSE = logFC / sqrt(LR)}, based on the likelihood ratio test statistic.
#'
#' When providing multiple contrasts, FDR correction can be applied globally using the
#' \code{global_fdr = TRUE} option.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param method A string: either \code{"edgeR"} or \code{"DESeq2"}.
#' @param group_var A string specifying the grouping variable in the sample metadata.
#'                 Not required if \code{contrast} is a named list for multiple variables.
#' @param contrast One of the following:
#'   \itemize{
#'     \item A character vector of two levels to compare (e.g., \code{c("Control", "Treated")}).
#'     \item A list of such character vectors for multiple contrasts within one grouping variable.
#'     \item A named list of such lists (e.g., \code{list(Treatment = list(c("A", "B"), c("A", "C")), Genotype = list(...))}).
#'   }
#' @param pseudocount A numeric value used to replace zero or negative counts (default = \code{1}).
#' @param significance_level A numeric value for the FDR threshold to determine significance (default = \code{0.05}).
#' @param output_csv_path Optional path to save results as CSV files. For multiple contrasts, each is saved separately.
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa (default = \code{"Genus"}).
#' @param palette A vector of colors for plotting significant and non-significant points (default = \code{c("#FFEB3B", "#073B4C")}).
#' @param global_fdr Logical. If \code{TRUE}, applies FDR correction across all contrasts (default = \code{FALSE}).
#' @return A list of result objects, or a single result if a single contrast is given. Each result includes:
#' \item{results}{A `data.frame` of differential abundance results with taxonomic annotations.}
#' \item{obj_significant}{A filtered `phyloseq` or `TreeSummarizedExperiment` object with significant taxa.}
#' \item{plot}{A `ggplot2` volcano plot.}
#' \item{bar_plot}{A `ggplot2` bar plot of log fold changes and standard errors.}
#' \item{bar_abundance_plot}{A `ggplot2` bar plot showing relative abundance of significant taxa across groups.}
#' @importFrom utils write.csv
#' @importFrom SummarizedExperiment mcols
#' @importFrom DESeq2 dispersions<-
#' @importFrom stats p.adjust relevel as.formula model.matrix
#' @importFrom grDevices col2rgb rgb
#' @importFrom dplyr rename mutate left_join filter group_by summarise distinct if_else
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_text geom_bar geom_linerange
#'   geom_vline coord_flip labs theme_minimal theme element_text scale_color_manual
#'   scale_fill_manual position_dodge element_blank element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom phyloseq transform_sample_counts tax_glom tax_table sample_data prune_taxa
#'   otu_table sample_names
#' @importFrom microbiome meta
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @importFrom limma makeContrasts
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq estimateSizeFactors estimateDispersionsGeneEst
#'   nbinomWaldTest results
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Single contrast (DESeq2)
#'   res_single <- perform_and_visualize_DA(
#'     obj = physeq_16SOTU,
#'     method = "DESeq2",
#'     group_var = "Diet",
#'     contrast = c("Insectivore", "Carnivore"),
#'     target_glom = "Genus"
#'   )
#'
#'   # Multiple contrasts for one group_var
#'   contrast_list <- list(
#'     c("Insectivore", "Carnivore"),
#'     c("Omnivore", "Herbivore")
#'   )
#'
#'   res_multi <- perform_and_visualize_DA(
#'     obj = physeq_16SOTU,
#'     method = "DESeq2",
#'     group_var = "Diet",
#'     contrast = contrast_list,
#'     global_fdr = TRUE
#'   )
#'
#'   # Named contrast list for multiple group_vars
#'   contrast_named <- list(
#'     Diet = list(
#'       c("Insectivore", "Carnivore"),
#'       c("Omnivore", "Carnivore")
#'     ),
#'     Animal.type = list(
#'       c("Frog", "Salamander")
#'     )
#'   )
#'
#'   res_multi_factor <- perform_and_visualize_DA(
#'     obj = physeq_16SOTU,
#'     method = "DESeq2",
#'     significance_level = 0.01,
#'     contrast = contrast_named,
#'     target_glom = "Genus",
#'     global_fdr = TRUE
#'   )
#'
#'   # Create a combined factor of Animal.type and Diet
#'   phyloseq::sample_data(physeq_16SOTU)$ComboGroup <- factor(interaction(
#'     phyloseq::sample_data(physeq_16SOTU)$Animal.type,
#'     phyloseq::sample_data(physeq_16SOTU)$Diet,
#'     drop = TRUE
#'   ))
#'
#'   # Define valid contrasts with sufficient sample sizes
#'   contrast_list <- list(
#'     c("Salamander.Insectivore", "Lizard.Insectivore"),
#'     c("Salamander.Carnivore", "Snake.Carnivore"),
#'     c("Salamander.Carnivore", "Frog.Carnivore")
#'   )
#'
#'   # Perform multi-contrast DA analysis on the ComboGroup
#'   res_combo <- perform_and_visualize_DA(
#'     obj = physeq_16SOTU,
#'     method = "DESeq2",
#'     group_var = "ComboGroup",
#'     contrast = contrast_list,
#'     target_glom = "Genus",
#'     global_fdr = TRUE
#'   )
#'
#'   # Example: view one of the plots
#'   res_combo$Frog.Insectivore_vs_Lizard.Insectivore$bar_plot
#' }
#'
#' @export
perform_and_visualize_DA <- function(obj,
                                     method,
                                     group_var,
                                     contrast,
                                     pseudocount = 1,
                                     significance_level = 0.05,
                                     output_csv_path = NULL,
                                     target_glom = "Genus",
                                     palette = c("#FFEB3B", "#073B4C"),
                                     global_fdr = TRUE) {
  # Detect if input is TSE and convert to phyloseq
  is_TSE <- inherits(obj, "TreeSummarizedExperiment")
  if (is_TSE) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  # Extract metadata
  metadata <- as.data.frame(microbiome::meta(obj))

  # ----------- NAMED MULTI-FACTOR CONTRASTS -----------------------
  if (is.list(contrast) && !is.null(names(contrast))) {
    results_all <- list()
    all_pvals <- c()
    all_keys <- c()

    for (gvar in names(contrast)) {
      contrast_list <- contrast[[gvar]]

      for (pair in contrast_list) {
        if (!is.character(pair) || length(pair) != 2) {
          stop("Error: each contrast must be a character vector of length 2.")
        }

        res <- perform_and_visualize_DA(
          obj = obj,
          method = method,
          group_var = gvar,
          contrast = pair,
          pseudocount = pseudocount,
          significance_level = significance_level,
          output_csv_path = NULL,
          target_glom = target_glom,
          palette = palette,
          global_fdr = FALSE
        )

        key <- paste0(gvar, "_", make.names(pair[1]), "_vs_", make.names(pair[2]))
        results_all[[key]] <- res

        if (!is.null(res$results)) {
          all_pvals <- c(all_pvals, res$results$pvalue)
          all_keys <- c(all_keys, rep(key, nrow(res$results)))
        }
      }
    }

    # Apply global FDR correction
    if (global_fdr && length(all_pvals) > 0) {
      global_adj <- p.adjust(all_pvals, method = "BH")
      index <- 1
      for (key in names(results_all)) {
        n_rows <- nrow(results_all[[key]]$results)
        results_all[[key]]$results$padj_global <- global_adj[index:(index + n_rows - 1)]
        results_all[[key]]$results$Significance_global <- ifelse(
          results_all[[key]]$results$padj_global < significance_level,
          "Significant", "Not Significant"
        )
        index <- index + n_rows
      }
    }

    return(results_all)
  }

  # ----- Validate contrast if not using named contrast -----
  if (is.null(group_var)) {
    stop("group_var must be specified if contrast is not a named list.")
  }


  # ----------- MULTI-CONTRAST SAME GROUP_VAR ---------------------
  if (is.list(contrast) && all(sapply(contrast, function(x) length(x) == 2))) {
    results_multi <- lapply(contrast, function(c_pair) {
      perform_and_visualize_DA(
        obj = obj,
        method = method,
        group_var = group_var,
        contrast = c_pair,
        pseudocount = pseudocount,
        significance_level = significance_level,
        output_csv_path = NULL,
        target_glom = target_glom,
        palette = palette,
        global_fdr = FALSE
      )
    })

    names(results_multi) <- sapply(contrast, function(c_pair) {
      paste0(make.names(c_pair[1]), "_vs_", make.names(c_pair[2]))
    })

    if (global_fdr) {
      all_pvals <- unlist(lapply(results_multi, function(x) x$results$pvalue))
      all_adj <- p.adjust(all_pvals, method = "BH")
      index <- 1
      for (nm in names(results_multi)) {
        n <- nrow(results_multi[[nm]]$results)
        results_multi[[nm]]$results$padj_global <- all_adj[index:(index + n - 1)]
        results_multi[[nm]]$results$Significance_global <- ifelse(
          results_multi[[nm]]$results$padj_global < significance_level,
          "Significant", "Not Significant"
        )
        index <- index + n
      }
    }

    if (!is.null(output_csv_path)) {
      dir.create(dirname(output_csv_path), recursive = TRUE, showWarnings = FALSE)
      for (nm in names(results_multi)) {
        utils::write.csv(results_multi[[nm]]$results,
          file.path(dirname(output_csv_path), paste0("DA_", nm, ".csv")),
          row.names = FALSE
        )
      }
    }

    return(results_multi)
  }
  # ----------- SINGLE CONTRAST STARTS BELOW ----------------------
  # Validate contrast format
  if (!is.character(contrast) || length(contrast) != 2) {
    stop("Error: contrast must be a character vector of length 2 at this stage.")
  }

  # ----------- SINGLE CONTRAST STARTS BELOW ----------------------

  # Clean group factor
  metadata[[group_var]] <- factor(metadata[[group_var]])
  levels(metadata[[group_var]]) <- make.names(gsub("-", "_", levels(metadata[[group_var]])))
  phyloseq::sample_data(obj) <- phyloseq::sample_data(metadata)

  # Validate contrast format (after multi-contrast logic is bypassed)
  if (!is.character(contrast) || length(contrast) != 2) {
    stop("Error: contrast must be a character vector of length 2 at this stage.")
  }

  # Filter counts and aggregate
  obj <- phyloseq::transform_sample_counts(obj, function(x) ifelse(x <= 0, pseudocount, x))
  obj <- phyloseq::tax_glom(obj, taxrank = target_glom)
  taxonomy_table <- as.data.frame(phyloseq::tax_table(obj))
  taxonomy_table$OTU <- rownames(taxonomy_table)

  # Process contrast
  contrast_fixed <- make.names(gsub("-", "_", contrast))
  available_levels <- levels(metadata[[group_var]])
  if (!all(contrast_fixed %in% available_levels)) {
    stop("Error: Contrast levels not found: ", paste(contrast_fixed, collapse = ", "))
  }

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
    res <- edgeR::topTags(lrt, n = Inf)$table |>
      dplyr::rename(pvalue = PValue) |>
      dplyr::mutate(
        FDR = p.adjust(pvalue, method = "BH"),
        padj = FDR,
        lfcSE = logFC / sqrt(LR),
        lfcSE = ifelse(is.nan(lfcSE) | is.infinite(lfcSE), NA, lfcSE),
        Significance = ifelse(FDR < significance_level, "Significant", "Not Significant"),
        group = ifelse(logFC > 0, contrast[2], contrast[1])
      )
    res$OTU <- rownames(res)
    res <- dplyr::left_join(res, taxonomy_table, by = "OTU")
    return(res)
  }

  perform_DESeq2 <- function(obj, group_var, contrast_fixed, significance_level, taxonomy_table) {
    otu <- round(as(phyloseq::otu_table(obj), "matrix"))
    metadata <- as.data.frame(microbiome::meta(obj))
    sample_ids <- phyloseq::sample_names(obj)
    metadata <- metadata[rownames(metadata) %in% sample_ids, , drop = FALSE]
    metadata <- metadata[match(sample_ids, rownames(metadata)), , drop = FALSE]
    metadata[[group_var]] <- stats::relevel(metadata[[group_var]], ref = contrast_fixed[1])

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, design = as.formula(paste("~", group_var)))

    dds <- tryCatch(
      DESeq2::DESeq(dds, sfType = "poscounts"),
      error = function(e) {
        message("DESeq2: Falling back to gene-wise dispersion estimates due to curve fitting failure.")
        dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- DESeq2::nbinomWaldTest(dds)
        return(dds)
      }
    )

    res <- as.data.frame(DESeq2::results(dds, contrast = c(group_var, contrast_fixed[2], contrast_fixed[1])))
    res <- res |>
      dplyr::rename(logFC = log2FoldChange, pvalue = pvalue) |>
      dplyr::mutate(
        FDR = p.adjust(pvalue, method = "BH"),
        padj = FDR,
        Significance = ifelse(padj < significance_level, "Significant", "Not Significant"),
        group = ifelse(logFC > 0, contrast[2], contrast[1])
      )

    res$OTU <- rownames(res)
    res <- dplyr::left_join(res, taxonomy_table, by = "OTU")
    return(res)
  }

  results <- if (method == "edgeR") {
    perform_edgeR(obj, group_var, contrast_fixed, significance_level, taxonomy_table)
  } else {
    perform_DESeq2(obj, group_var, contrast_fixed, significance_level, taxonomy_table)
  }

  results <- results |>
    dplyr::filter(!is.na(Significance)) |>
    dplyr::mutate(group = factor(group))
  sig_otus <- results$OTU[results$Significance == "Significant"]
  if (length(sig_otus) > 0) {
    obj_significant <- phyloseq::prune_taxa(sig_otus, obj)
  } else {
    obj_significant <- NULL
    message("No significant taxa found for contrast: ", paste(contrast, collapse = " vs "))
  }
  p <- if (nrow(results) > 0) {
    ggplot2::ggplot(results, ggplot2::aes(
      x = logFC, y = -log10(pvalue),
      color = Significance,
      shape = group,
      size = -log10(FDR)
    )) +
      ggplot2::geom_point(alpha = 0.8) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = ifelse(FDR < significance_level & -log10(pvalue) > 5, OTU, "")),
        size = 4, fontface = "plain", max.overlaps = 15, segment.color = "grey40"
      ) +
      ggplot2::geom_hline(yintercept = -log10(significance_level), linetype = "dashed", color = "#FF5733", linewidth = 1) +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::labs(
        x = "Log2 Fold Change",
        y = "-log10 (P-value)",
        title = paste("Volcano Plot -", group_var),
        subtitle = paste("FDR threshold:", significance_level)
      )
  } else {
    NULL
  }

  df_filtered <- results |>
    dplyr::group_by(Genus, group) |>
    dplyr::summarise(
      logFC = mean(logFC, na.rm = TRUE),
      lfcSE = mean(lfcSE, na.rm = TRUE),
      padj = ifelse(all(is.na(padj)), NA, min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::filter(!is.na(Genus) & !is.na(padj) & padj < significance_level)
  df_filtered <- df_filtered |> dplyr::distinct(Genus, .keep_all = TRUE)
  df_filtered$Genus <- factor(df_filtered$Genus, levels = unique(df_filtered$Genus[order(df_filtered$logFC)]))
  df_filtered <- df_filtered |> dplyr::mutate(
    LFC_Direction = ifelse(logFC < 0, "Negative LFC", "Positive LFC"),
    group_label = paste0(group, " (", LFC_Direction, ")")
  )

  group_levels <- unique(as.character(results$group))
  fixed_group_colors <- c("#9183E6", "#33FFD1")[seq_along(group_levels)]
  names(fixed_group_colors) <- group_levels

  bar_colors <- NULL
  linerange_colors <- NULL

  if (nrow(df_filtered) > 0) {
    bar_colors <- setNames(fixed_group_colors[df_filtered$group], df_filtered$group_label)
    bar_colors <- bar_colors[!duplicated(names(bar_colors))]

    darken_hex <- function(hex_vec, factor = 0.55) {
      hex_vec <- as.character(hex_vec)
      rgb_mat <- grDevices::col2rgb(hex_vec) / 255
      rgb_dark <- apply(rgb_mat, 2, function(col) pmax(0, col * factor))
      grDevices::rgb(rgb_dark[1, ], rgb_dark[2, ], rgb_dark[3, ])
    }

    linerange_colors <- setNames(darken_hex(bar_colors), names(bar_colors))
  }

  bar_plot <- NULL
  bar_abundance_plot <- NULL

  if (nrow(df_filtered) > 0) {
    bar_plot <- ggplot2::ggplot(
      df_filtered,
      ggplot2::aes(
        x = Genus,
        y = logFC,
        fill = group_label,
        color = group_label
      )
    ) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(0.7), alpha = 0.9, width = 0.8) +
      ggplot2::geom_point(size = 0.6, position = ggplot2::position_dodge(0.7), shape = 21, stroke = 0.6) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = logFC - lfcSE, ymax = logFC + lfcSE),
        size = 0.9
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        subtitle = paste("Significant taxa with padj <", significance_level),
        caption = "Bars = log2 fold change; Lines = +/- 1 SE",
        x = "",
        y = "Log Fold Change",
        fill = "Group & LFC Direction",
        color = "Group & LFC Direction"
      ) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5, color = "gray30"),
        axis.text.y = ggplot2::element_text(face = "italic", color = "black", size = 12),
        axis.text.x = ggplot2::element_text(face = "bold", size = 12),
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "top",
        legend.title = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_manual(values = bar_colors) +
      ggplot2::scale_color_manual(values = linerange_colors)

    bar_abundance_plot <- plotbar_abundance(
      physeq = obj_significant,
      tax_level = target_glom,
      normalize = TRUE,
      treatment_variable = group_var,
      abundance_type = "absolute",
      palette = DspikeIn::color_palette$mix_MG
    )
  }


  if (is_TSE) {
    obj_significant <- convert_phyloseq_to_tse(obj_significant)
  }

  return(list(
    results = results,
    obj_significant = obj_significant,
    plot = p,
    bar_plot = bar_plot,
    bar_abundance_plot = bar_abundance_plot
  ))
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
# results_edgeR$bar_plot

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
# results_DESeq2$bar_plot
