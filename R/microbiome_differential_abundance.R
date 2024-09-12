#' Install and Load Required Packages
#' This function installs and loads required packages.
#' 
#' @param packages A vector of package names to be installed and loaded.
#' @examples
#' required_packages <- c("phyloseq", "DESeq2", "edgeR", "BiocManager", "BiocGenerics", "ggplot2", "dplyr", "DT")
#' install_and_load(required_packages)
#' @export
#' @importFrom utils install.packages
#' @importFrom BiocManager install
install_and_load <- function(packages) {
  for (package in packages) {
    # Check if the package is already installed
    if (!requireNamespace(package, quietly = TRUE)) {
      suppressMessages({
        suppressWarnings({
          # Check if it is a Bioconductor package
          if (package %in% c("phyloseq", "DESeq2", "edgeR", "BiocManager", "BiocGenerics")) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              utils::install.packages("BiocManager")
            }
            BiocManager::install(package, suppressUpdates = TRUE, ask = FALSE)
          } else {
            utils::install.packages(package)
          }
        })
      })
    }
    # Load the package
    suppressMessages(suppressWarnings(library(package, character.only = TRUE)))
  }
}

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
#' @param ps A phyloseq object.
#' @param pseudocount A numeric value to add to avoid zero counts.
#' @return A phyloseq object with filtered and adjusted OTU table.
#' @examples
#' ps <- remove_zero_negative_count_samples(ps)
#' @export
remove_zero_negative_count_samples <- function(ps, pseudocount = 1) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  
  zero_negative_count_samples <- phyloseq::sample_sums(ps) <= 0
  na_count_samples <- apply(otu, 2, function(x) any(is.na(x)))
  samples_to_remove <- zero_negative_count_samples | na_count_samples
  if (any(samples_to_remove)) {
    cat("Removing", sum(samples_to_remove), "samples with zero, negative counts, or NA values.\n")
    ps <- phyloseq::prune_samples(!samples_to_remove, ps)
    otu <- as(phyloseq::otu_table(ps), "matrix")
  }
  
  zero_rows <- rowSums(otu) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples.\n")
    otu <- otu[!zero_rows, ]
    ps <- phyloseq::prune_taxa(!zero_rows, ps)
  }
  
  otu <- otu + pseudocount
  otu <- round(otu)
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  return(ps)
}

#' Convert Categorical Columns to Factors in Sample Data
#' @param ps A phyloseq object.
#' @return A phyloseq object with updated sample data.
#' @examples
#' ps <- convert_categorical_to_factors(ps)
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

#' Relativized Filtered Taxa
#' This function filters taxa from a phyloseq object based on custom thresholds.
#' @param physeq A phyloseq object containing the microbial data.
#' @param threshold_percentage A numeric value specifying the minimum percentage of samples in which a taxon must be present.
#' @param threshold_mean_abundance A numeric value specifying the minimum mean abundance.
#' @param threshold_count A numeric value specifying the minimum count of a taxon in a sample.
#' @param threshold_relative_abundance A numeric value specifying the minimum relative abundance.
#' @return A phyloseq object containing only the taxa that meet the specified thresholds.
#' @examples
#' FT <- relativized_filtered_taxa(spiked_16S, threshold_percentage = 0.6, threshold_mean_abundance = 0.0005, threshold_count = 5)
#' @export
relativized_filtered_taxa <- function(physeq, threshold_percentage = 0.5, threshold_mean_abundance = 0.001, threshold_count = 10, threshold_relative_abundance = NULL) {
  nsamples <- phyloseq::nsamples(physeq)
  sample_sum <- phyloseq::sample_sums(physeq)
  
  filter_function <- function(x) {
    (sum(x > threshold_count) > nsamples * threshold_percentage) | 
      ((sum(x > threshold_count) > (nsamples * 0.1)) & (mean(x / sample_sum) > threshold_mean_abundance) & (max(x / sample_sum) > threshold_relative_abundance))
  }
  
  two_way_filtered <- phyloseq::filter_taxa(physeq, filter_function, prune = TRUE)
  return(two_way_filtered)
}

#' #Glom Taxa at Specific Level
#' #This function gloms taxa at a specified taxonomic rank. If no rank is provided, it defaults to "Genus".
#' @param physeq A phyloseq object.
#' @param rank A character string specifying the taxonomic rank to glom. Default is "Genus".
#' @return A phyloseq object with taxa glommed at the specified rank.
#' @examples
#' ps_glommed <- glom_taxa_at_rank(ps) # Defaults to "Genus"
#' ps_glommed <- glom_taxa_at_rank(ps, "Family") # Gloms at "Family"
#' @export
glom_taxa_at_rank <- function(physeq, rank = "Genus") {
  if (is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
    stop("The tax_glom() function requires that physeq contain a taxonomyTable")
  }
  physeq <- phyloseq::tax_glom(physeq, taxrank = rank)
  
  if (phyloseq::ntaxa(physeq) == 0) stop("No taxa remain after glomming.")
  
  return(physeq)
}

#' Perform Differential Abundance Analysis with DESeq2
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#'
#' @param ps A phyloseq object containing the microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage. Default is 0.0001.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance. Default is 0.00001.
#' @param threshold_count A numeric value for filtering count. Default is 1.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance. Default is 0.00001.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 3.
#' @param facet_variable A string specifying the variable to facet the barplots. Default is "Phylum".
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa (e.g., "Genus", "Family"). Default is "Genus".
#' @return A list containing the final results, the phyloseq object with significant OTUs, and ggplot objects for the visualizations.
#' @export
perform_and_visualize_differential_abundance <- function(ps, group_var, threshold_percentage = 0.0001, 
                                                         threshold_mean_abundance = 0.00001, threshold_count = 1, 
                                                         threshold_relative_abundance = 0.00001, significance_level = 0.05, 
                                                         point_size = 3, facet_variable = "Phylum", target_glom = "Genus") {
  
  # Step 1: Prune samples and taxa with zero counts
  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 0, ps)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  
  # Step 2: Check if the phyloseq object has a taxonomy table before glomming
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    ps <- glom_taxa_at_rank(ps, target_glom)
  } else {
    warning("No taxonomy table found. Skipping tax glomming step.")
  }
  
  # Step 3: Perform DESeq2 differential abundance analysis
  results <- perform_DESeq2(ps, group_var, threshold_percentage, threshold_mean_abundance, 
                            threshold_count, threshold_relative_abundance, significance_level)
  
  # Rename DESeq2 log2FoldChange to logFC for consistency
  results <- results %>%
    dplyr::rename(logFC = log2FoldChange)
  
  # Step 4: Create Volcano Plot
  volcano_plot <- ggplot(results, aes(x = logFC, y = pvalue)) +
    geom_point(aes(color = diff_abn), size = point_size) +
    scale_y_continuous(trans = "log10", labels = scales::label_scientific()) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "P-value", title = "Volcano Plot") +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "bottom")
  
  # Step 5: Create Barplots for Relative and Absolute Abundances
  barplot_rel <- taxa_barplot(ps, target_glom = target_glom, treatment_variable = group_var, 
                              abundance_type = "relative", x_angle = 90, fill_variable = "Genus", 
                              facet_variable = facet_variable, top_n_taxa = 20, palette = MG())$barplot
  
  barplot_abs <- taxa_barplot(ps, target_glom = target_glom, treatment_variable = group_var, 
                              abundance_type = "absolute", x_angle = 90, fill_variable = "Genus", 
                              facet_variable = facet_variable, top_n_taxa = 20, palette = MG())$barplot
  
  # Return results and plots
  return(list(plot = volcano_plot, barplot_rel = barplot_rel, barplot_abs = barplot_abs, results = results))
}

#' Glom Taxa at Specific Level
#' This function gloms taxa at a specified taxonomic rank.
#' @param physeq A phyloseq object.
#' @param rank A character string specifying the taxonomic rank to glom. Default is "Genus".
#' @return A phyloseq object with taxa glommed at the specified rank.
#' @export
glom_taxa_at_rank <- function(physeq, rank = "Genus") {
  if (is.null(phyloseq::tax_table(physeq, errorIfNULL = FALSE))) {
    stop("The tax_glom() function requires that physeq contain a taxonomyTable")
  }
  physeq <- phyloseq::tax_glom(physeq, taxrank = rank)
  
  if (phyloseq::ntaxa(physeq) == 0) stop("No taxa remain after glomming.")
  
  return(physeq)
}
#' Perform Differential Abundance Analysis with DESeq2
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param rank A string specifying the taxonomic rank to glom (e.g., "Genus", "Family"). Default is "Genus".
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
perform_DESeq2 <- function(ps, group_var, rank = "Genus", threshold_percentage, threshold_mean_abundance, 
                           threshold_count, threshold_relative_abundance, significance_level) {
  
  # Step 1: Glom taxa at the specified rank (e.g., Genus)
  ps_glommed <- glom_taxa_at_rank(ps, rank = rank)
  
  # Step 2: Prune samples and taxa with zero counts after glomming
  ps_glommed <- phyloseq::prune_samples(phyloseq::sample_sums(ps_glommed) > 0, ps_glommed)
  ps_glommed <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_glommed) > 0, ps_glommed)
  
  # Step 3: Ensure sample names are aligned
  if (!all(colnames(phyloseq::otu_table(ps_glommed)) == rownames(phyloseq::sample_data(ps_glommed)))) {
    stop("Sample names in OTU table and sample data do not match.")
  }
  
  # Convert the phyloseq object to a DESeq2-compatible matrix
  otu <- as(phyloseq::otu_table(ps_glommed), "matrix")
  metadata <- as(phyloseq::sample_data(ps_glommed), "data.frame")
  
  # Check if the group variable exists in metadata
  if (!group_var %in% colnames(metadata)) {
    stop("The specified group variable does not exist in the sample data.")
  }
  
  # Filter taxa using custom thresholds
  ps_filtered <- tryCatch({
    relativized_filtered_taxa(ps_glommed, threshold_percentage, threshold_mean_abundance, 
                              threshold_count, threshold_relative_abundance)
  }, error = function(e) {
    warning("Filtering resulted in no taxa remaining: ", e$message, ". Skipping filtering step.")
    return(ps_glommed)
  })
  
  if (phyloseq::nsamples(ps_filtered) == 0 || phyloseq::ntaxa(ps_filtered) == 0) {
    stop("No taxa remain after filtering. Adjust the filtering thresholds.")
  }
  
  # Glom taxa at the rank specified
  ps_filtered_glom <- glom_taxa_at_rank(ps_filtered, rank)
  
  # Round the OTU table
  otu <- as(phyloseq::otu_table(ps_filtered_glom), "matrix")
  otu <- round(otu)
  phyloseq::otu_table(ps_filtered_glom) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  
  # Create DESeq2 dataset and run DESeq
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = otu, colData = metadata, 
                                        design = as.formula(paste("~", group_var)))
  dds <- DESeq2::DESeq(dds)
  
  # Extract results and arrange by p-value
  res <- DESeq2::results(dds)
  res <- as.data.frame(res)
  res$OTU <- rownames(res)
  
  # Adjust for multiple testing using FDR and mark significant results
  res <- res %>%
    dplyr::arrange(pvalue) %>%
    dplyr::mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    dplyr::mutate(diff_abn = FDR < significance_level)
  
  return(res)
}
#' Perform Differential Abundance Analysis with edgeR
#' @param ps A phyloseq object.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param threshold_percentage A numeric value for filtering threshold percentage.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance.
#' @param threshold_count A numeric value for filtering count.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs.
#' @return A data frame containing the differential abundance results.
#' @export
#' @importFrom edgeR DGEList estimateDisp glmFit glmLRT topTags
perform_edgeR <- function(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level) {
  # Extract OTU table and metadata from phyloseq object
  otu <- as(otu_table(ps), "matrix")
  metadata <- as(sample_data(ps), "data.frame")
  
  # Check if the grouping variable exists in the sample data
  if (!group_var %in% colnames(metadata)) {
    stop("The specified group variable does not exist in the sample data.")
  }
  
  # Filter taxa based on the thresholds
  ps_filtered <- tryCatch({
    relativized_filtered_taxa(ps, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance)
  }, error = function(e) {
    warning("Filtering resulted in no taxa remaining: ", e$message, ". Skipping filtering step.")
    return(ps)
  })
  
  # Check if there are any remaining taxa or samples after filtering
  if (nsamples(ps_filtered) == 0 || ntaxa(ps_filtered) == 0) {
    stop("No taxa remain after filtering. Adjust the filtering thresholds.")
  }
  
  # Glom taxa at the genus rank
  ps_filtered_glom <- glom_taxa_at_rank(ps_filtered, "Genus")
  
  # Normalize the OTU table
  otu_normalized <- as(otu_table(ps_filtered_glom), "matrix")
  group <- metadata[[group_var]]
  
  # Create DGEList object for edgeR analysis
  dge <- edgeR::DGEList(counts = otu_normalized, group = group)
  
  # Design matrix
  design <- model.matrix(~ group)
  
  # Estimate dispersions
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit the model and perform likelihood ratio test
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit)
  
  # Extract top tags (significant results)
  top_tags <- edgeR::topTags(lrt, n = Inf)
  results <- top_tags$table
  
  # Add OTU names and compute FDR
  results$OTU <- rownames(results)
  results <- results %>%
    arrange(PValue) %>%
    mutate(FDR = p.adjust(PValue, method = "BH")) %>%
    mutate(diff_abn = FDR < significance_level)
  
  return(results)
}

#' Perform and Visualize Differential Abundance Analysis
#' This function normalizes and filters the data, performs differential abundance analysis, 
#' extracts significant OTUs, rebuilds the phyloseq object, merges results with metadata, and visualizes the results.
#' 
#' @param ps A phyloseq object containing the microbiome data.
#' @param group_var A string specifying the grouping variable in sample data.
#' @param method A string specifying the differential abundance method ("edgeR", "DESeq2"). Default is "edgeR".
#' @param threshold_percentage A numeric value for filtering threshold percentage. Default is 0.0001.
#' @param threshold_mean_abundance A numeric value for filtering mean abundance. Default is 0.00001.
#' @param threshold_count A numeric value for filtering count. Default is 1.
#' @param threshold_relative_abundance A numeric value for filtering relative abundance. Default is 0.00001.
#' @param significance_level A numeric value specifying the significance level for filtering significant OTUs. Default is 0.05.
#' @param point_size A numeric value specifying the size of points in the volcano plot. Default is 3.
#' @param facet_variable A string specifying the variable to facet the barplots. Default is "Phylum".
#' @param target_glom A string specifying the taxonomic rank to aggregate taxa (e.g., "Genus", "Family"). Default is "Genus".
#' @return A list containing the final results, the phyloseq object with significant OTUs, and ggplot objects for the visualizations.
#' 
#' @examples
#' # Install and load the required packages first
#' required_packages <- c("phyloseq", "DESeq2", "edgeR", "BiocManager", "BiocGenerics", "ggplot2", "dplyr", "DT")
#' install_and_load(required_packages)
#' 
#' # Step 1: Clean and preprocess the data
#' ps <- remove_zero_negative_count_samples(ps)
#' ps <- convert_categorical_to_factors(ps)
#' 
#' # edgeR Example with default glomming at the "Genus" rank:
#' results_edgeR <- perform_and_visualize_differential_abundance(ps, "Animal.ecomode", method = "edgeR", 
#'                                                               threshold_percentage = 0.001, threshold_mean_abundance = 0.001,
#'                                                               threshold_count = 5, threshold_relative_abundance = 0.001,
#'                                                               significance_level = 0.05, point_size = 3, target_glom = "Genus")
#' # Print the volcano plot and barplots for edgeR
#' print(results_edgeR$plot)
#' print(results_edgeR$barplot_rel)
#' print(results_edgeR$barplot_abs)
#' 
#' # DESeq2 Example with custom glomming at the "Family" rank and faceting by Class:
#' results_DESeq2 <- perform_and_visualize_differential_abundance(ps, "Host.genus", method = "DESeq2",
#'                                                                threshold_percentage = 0.001, threshold_mean_abundance = 0.001,
#'                                                                threshold_count = 5, threshold_relative_abundance = 0.001,
#'                                                                significance_level = 0.05, point_size = 3, target_glom = "Family", 
#'                                                                facet_variable = "Class")
#' # Print the volcano plot and barplots for DESeq2
#' print(results_DESeq2$plot)
#' print(results_DESeq2$barplot_rel)
#' print(results_DESeq2$barplot_abs)
#' 
#' @export
perform_and_visualize_differential_abundance <- function(ps, group_var, method = "edgeR",
                                                         threshold_percentage = 0.0001, threshold_mean_abundance = 0.00001, 
                                                         threshold_count = 1, threshold_relative_abundance = 0.00001,
                                                         significance_level = 0.05, point_size = 3, facet_variable = "Phylum",
                                                         target_glom = "Genus") {
  
  # Step 1: Check if the grouping variable has at least two levels
  group_levels <- levels(as.factor(phyloseq::sample_data(ps)[[group_var]]))
  if (length(group_levels) < 2) {
    stop(paste("The grouping variable '", group_var, "' must have at least two levels. Currently, it has the following levels: ", paste(group_levels, collapse = ", "), sep = ""))
  }
  
  # Step 2: Prune samples and taxa with zero counts
  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 0, ps)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  
  # Step 3: Check if the physeq object has a taxonomy table before glomming
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    ps <- glom_taxa_at_rank(ps, target_glom)  # Proceed with tax glom at user-specified rank
  } else {
    warning("No taxonomy table found. Skipping tax glomming step.")
  }
  
  # Step 4: Perform differential abundance analysis based on the chosen method
  if (method == "edgeR") {
    results <- perform_edgeR(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level)
  } else if (method == "DESeq2") {
    results <- perform_DESeq2(ps, group_var, threshold_percentage, threshold_mean_abundance, threshold_count, threshold_relative_abundance, significance_level)
  } else {
    stop("Invalid method specified. Use 'edgeR' or 'DESeq2'.")
  }
  
  # Adjust for DESeq2 logFC column
  if (method == "DESeq2") {
    results <- results %>%
      dplyr::rename(logFC = log2FoldChange)  # DESeq2 uses log2FoldChange instead of logFC
  }
  
  # Step 5: Generate visualizations based on the analysis results
  p <- visualize_differential_abundance(results, group_var, point_size = point_size, palette = MG())
  
  bp_rel <- taxa_barplot(ps, target_glom = target_glom, treatment_variable = group_var, 
                         abundance_type = "relative", x_angle = 90, fill_variable = "Genus", facet_variable = facet_variable, top_n_taxa = 20, palette = MG())
  
  bp_abs <- taxa_barplot(ps, target_glom = target_glom, treatment_variable = group_var, 
                         abundance_type = "absolute", x_angle = 90, fill_variable = "Genus", facet_variable = facet_variable, top_n_taxa = 20, palette = MG())
  
  return(list(plot = p, barplot_rel = bp_rel$barplot, barplot_abs = bp_abs$barplot, results = results))
}


# Example usage:
# Assuming 'ps' is a phyloseq object with group variables
# Remove samples with zero, negative, or NA counts and add a pseudocount
# ps <- remove_zero_negative_count_samples(ps)
# 
# # Step 2: Ensure that categorical variables in the sample data are factors
# ps <- convert_categorical_to_factors(ps)

# Step 3: Perform differential abundance analysis and visualize results
# edgeR/DESeq2 Example with default faceting by "Phylum":
# Diet<-as.factor(ps@sam_data$Diet)
# edgeR Example with default glomming at the "Genus" rank
# results_edgeR <- perform_and_visualize_differential_abundance(ps, 
#                                                               group_var = "Animal.ecomode", 
#                                                               method = "edgeR", 
#                                                               threshold_percentage = 0.001, 
#                                                               threshold_mean_abundance = 0.001,
#                                                               threshold_count = 5, 
#                                                               threshold_relative_abundance = 0.001,
#                                                               significance_level = 0.05, 
#                                                               point_size = 3, facet_variable = "Animal.ecomode",
#                                                               target_glom = "Genus")
# 
# # Print results from edgeR
# print(results_edgeR$plot)
# print(results_edgeR$barplot_rel)
# print(results_edgeR$barplot_abs)
# 
# # DESeq2 Example with custom glomming at the "Family" rank and faceting by Class
# results_DESeq2 <- perform_and_visualize_differential_abundance(ps,
#                                                                group_var = "Diet",
#                                                                threshold_percentage = 0.001,
#                                                                threshold_mean_abundance = 0.001,
#                                                                threshold_count = 5,
#                                                                threshold_relative_abundance = 0.001,
#                                                                significance_level = 0.05,
#                                                                point_size = 3,
#                                                                facet_variable = "Animal.ecomode")
# 
# # Print the volcano plot and barplots
# print(results_DESeq2$plot)
# print(results_DESeq2$barplot_rel)
# print(results_DESeq2$barplot_abs)


