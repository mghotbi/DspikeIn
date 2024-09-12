# For Differential Abundance Analysis: Consider DESeq, TMM, or CSS normalization, as these methods are designed to handle compositional biases and library size differences.
# For Compositional Data: CLR normalization is a good choice, as it transforms the data to account for its compositional nature.
# For Simplicity and Ease of Use: TC, UQ, or Median normalization are quick and straightforward but may not be as robust.
# Install and Load Required Packages 
#
# Differential Abundance Analysis:
# DESeq, TMM, or CSS normalization methods are designed to handle compositional biases 
# and library size differences effectively. 

# For Compositional Data:
# CLR normalization is ideal as it accounts for the compositional nature of the data 
# by transforming it into log-ratio format.
# For Simplicity and Ease of Use:
# TC, UQ, or Median normalization methods are quick and easy but might not be as robust.

#' -----------------------------------------------------------
#' @importFrom phyloseq otu_table sample_data sample_sums prune_samples prune_taxa taxa_are_rows tax_table sample_names
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq sizeFactors counts
#' @importFrom BiocManager install
#' @importFrom stats model.matrix median
#' 
# ------------------------------------------------------------------------
#'
#' @param packages A character vector of package names to install and load.
install_and_load <- function(packages) {
  for (package in packages) {
    suppressMessages({
      if (!requireNamespace(package, quietly = TRUE)) {
        if (package %in% c("phyloseq", "DESeq2", "edgeR", "EDASeq", "BiocManager", "BiocGenerics")) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            utils::install.packages("BiocManager", quiet = TRUE)
          }
          BiocManager::install(package, suppressUpdates = TRUE, ask = FALSE, quiet = TRUE)
        } else {
          utils::install.packages(package, quiet = TRUE)
        }
      }
      suppressWarnings(suppressMessages(library(package, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    })
  }
}

# Example usage
required_packages <- c("phyloseq", "DESeq2", "edgeR", "EDASeq", "BiocManager", "BiocGenerics")
install_and_load(required_packages)

# ------------------------------------------------------------------------

#' Calculate Geometric Mean
#'
#' @param x A numeric vector.
#' @param na.rm Logical, should missing values (NAs) be removed?
#' @return Geometric mean of x.
gm_mean <- function(x, na.rm = TRUE) {
  valid_x <- x[x > 0 & !is.na(x)]
  if (length(valid_x) == 0) return(NA)
  exp(sum(log(valid_x), na.rm = na.rm) / length(valid_x))
}

# ------------------------------------------------------------------------

#' Set Normalization Factors in the Sample Data of the Phyloseq Object
#'
#' @param ps A phyloseq object.
#' @param scaling.factor A vector of normalization factors.
#' @return A phyloseq object with updated sample data.
#' @export
set_nf <- function(ps, scaling.factor) {
  if (!inherits(ps, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  
  if (length(scaling.factor) != phyloseq::nsamples(ps)) {
    stop("Length of scaling.factor must match the number of samples in the phyloseq object.")
  }
  
  phyloseq::sample_data(ps)$norm_factors <- scaling.factor
  return(ps)
}

# -----------------------------------------------------------
#' Remove Samples with Zero, Negative Counts, or NA Values and Add Pseudocount
#' 
#' @param ps A phyloseq object.
#' @param pseudocount A numeric value to add to avoid zero counts.
#' @return A phyloseq object with filtered and adjusted OTU table.
remove_zero_negative_count_samples <- function(ps, pseudocount = 1e-6) {
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

# -----------------------------------------------------------
#' Convert Categorical Columns to Factors in Sample Data
#' 
#' @param ps A phyloseq object.
#' @return A phyloseq object with updated sample data.
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

# -----------------------------------------------------------
#' Create a List from a Phyloseq Object
#' 
#' @param physeq A phyloseq object.
#' @return A list containing the DGE list and updated phyloseq object.
create_list <- function(physeq) {
  if (!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  counts <- as(phyloseq::otu_table(physeq), "matrix")
  sample_data_df <- as(phyloseq::sample_data(physeq), "data.frame")
  taxonomy <- phyloseq::tax_table(physeq)
  lib.size <- colSums(counts)
  norm.factors <- rep(1, ncol(counts))
  sample_data_df$lib.size <- lib.size
  sample_data_df$norm.factors <- norm.factors
  dge_base <- list(
    counts = counts,
    samples = sample_data_df,
    genes = NULL,
    group = sample_data_df$group,  
    lib.size = lib.size,
    norm.factors = norm.factors
  )
  updated_physeq <- physeq
  phyloseq::sample_data(updated_physeq) <- phyloseq::sample_data(sample_data_df)
  return(list(dge_list = dge_base, phyloseq_obj = updated_physeq))
}

# -----------------------------------------------------------
#' Apply the Selected Normalization Method to the Phyloseq Object
#' 
#' @param ps A phyloseq object.
#' @param method A character string specifying the normalization method ("TC", "UQ", "med", "DESeq", "Poisson", "QN", "TMM", "clr", "rar", "css", "tss", "rle").
#' @param groups A column name of group labels from sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @export
normalization_set <- function(ps, method, groups = NULL) {
  if (phyloseq::nsamples(ps) == 0) stop("The phyloseq object contains no samples.")
  if (length(phyloseq::sample_names(ps)) == 0) stop("Sample names are missing in the phyloseq object.")
  
  ps <- remove_zero_negative_count_samples(ps)
  ps <- convert_categorical_to_factors(ps)
  
  # Ensure groups is a single column name from sample data
  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) != 1 || !groups %in% colnames(phyloseq::sample_data(ps))) {
      stop("'groups' should be a single column name from sample data.")
    }
  }
  
  result <- switch(method,
                         "TC" = norm.TC(ps, groups),
                         "UQ" = norm.UQ(ps, groups),
                         "med" = norm.med(ps, groups),
                         "DESeq" = norm.DESeq(ps, groups),
                         "Poisson" = norm.Poisson(ps, groups),
                         "QN" = norm.QN(ps),
                         "TMM" = norm.TMM(ps, groups),
                         "clr" = norm.clr(ps),
                         "rar" = norm.rar(ps),
                         "css" = norm.css(ps),
                         "tss" = norm.tss(ps),
                         "rle" = norm.rle(ps),
                         stop("Invalid normalization method"))
  
  dat.normed <- result$dat.normed
  scaling.factor <- result$scaling.factor
  
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' TC Normalization (Total Count Scaling)
#' 
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.TC <- function(ps, groups) {
  ps <- create_list(ps)$phyloseq_obj
  dat.DGE <- create_list(ps)$dge_list
  scaling.factor <- dat.DGE$samples$lib.size / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' UQ Normalization (Upper Quartile)
#' 
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.UQ <- function(ps, groups) {
  # Create a custom list
  physeq_list <- create_list(ps)
  ps <- physeq_list$phyloseq_obj
  dat.DGE <- physeq_list$dge_list
  
  # Calculate the upper quartile factor
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  
  # Handle cases where q.factor might be zero
  if (any(q.factor == 0)) stop("One or more upper quartile factors are zero.")
  
  scaling.factor <- q.factor / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  
  # Update the OTU table in the phyloseq object
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  
  # Set normalization factors
  ps <- set_nf(ps, scaling.factor)
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' Median Normalization
#' 
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.med <- function(ps, groups) {
  ps <- create_list(ps)$phyloseq_obj
  dat.DGE <- create_list(ps)$dge_list
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' DESeq Normalization
#' 
# -----------------------------------------------------------
#' DESeq Normalization with Pseudocount and Integer Conversion
#' 
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @param pseudocount A numeric value added to avoid zeros in the dataset.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.DESeq <- function(ps, groups, pseudocount = 1) {
  ps <- remove_zero_negative_count_samples(ps)  # Clean the data
  
  # Convert OTU table to a matrix and add pseudocount
  raw <- as(phyloseq::otu_table(ps), "matrix") + pseudocount
  
  # Round counts to integers
  raw <- round(raw)
  
  # Check if the groups vector has only one unique value
  unique_groups <- unique(groups)
  if (length(unique_groups) == 1) {
    # All samples have the same condition, use design ~ 1
    design <- stats::model.matrix(~ 1)
  } else {
    design <- stats::model.matrix(~ groups)
  }
  
  # Create condition data frame for DESeq2
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) <- colnames(raw)
  
  # Create DESeq2 dataset (counts must be integers)
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = condition, design = design)
  
  # Run DESeq normalization
  if (length(unique_groups) == 1) {
    dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  } else {
    dat.DGE <- DESeq2::DESeq(dat.DGE, fitType = "local")
  }
  
  # Get scaling factors and normalized data
  scaling.factor <- DESeq2::sizeFactors(dat.DGE)
  dat.normed <- DESeq2::counts(dat.DGE, normalized = TRUE)
  
  # Update phyloseq object with normalized data
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' Quantile Normalization (QN) for phyloseq object 
#' 
#' @param ps A phyloseq object.
#' @param filter Logical, whether to filter low counts.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.QN <- function(ps, filter = FALSE) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  
  if (!is.numeric(otu)) {
    stop("OTU table must contain numeric values.")
  }
  
  taxa_names_original <- phyloseq::taxa_names(ps)
  tax_table_original <- phyloseq::tax_table(ps)
  
  if (filter) {
    otu <- log2(otu + 1)
    otu <- otu[rowMeans(otu) > 2, ]
  } else {
    otu <- log2(otu + 1)
  }
  
  rank_mean <- apply(otu, 2, rank)
  sorted <- apply(otu, 2, sort)
  mean_values <- rowMeans(sorted)
  normalized <- apply(rank_mean, 2, function(r) mean_values[round(r)])
  normalized <- 2^normalized - 1
  rownames(normalized) <- taxa_names_original
  phyloseq::otu_table(ps) <- phyloseq::otu_table(normalized, taxa_are_rows = TRUE)
  
  # for QN we don't need scaling factors
  scaling.factor <- NULL
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------

#' Poisson Normalization and Differential Abundance Function
#' 
#' @param ps A phyloseq object or matrix of raw counts.
#' @param group_var A string specifying the grouping variable in sample data (if phyloseq object).
#' @param pseudocount A numeric value added to avoid division by zero.
#' @return A list containing the normalized data, scaling factor, and differential abundance results.
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
norm.Poisson <- function(ps, group_var = NULL, pseudocount = 1e-6) {
  ps <- remove_zero_negative_count_samples(ps)
  
  prepare_data <- function(ps, group_var) {
    if (inherits(ps, "phyloseq")) {
      raw_otu <- as(phyloseq::otu_table(ps), "matrix")
      sample_data_df <- as(phyloseq::sample_data(ps), "data.frame")
      
      if (!is.null(group_var)) {
        y <- as.numeric(as.factor(sample_data_df[[group_var]]))
      } else {
        y <- rep(1, ncol(raw_otu))
      }
    } else if (is.matrix(ps)) {
      raw_otu <- ps
      y <- rep(1, ncol(raw_otu))
    } else {
      stop("Input ps must be a phyloseq object or a matrix.")
    }
    list(raw_otu = raw_otu, sample_data_df = sample_data_df, y = y)
  }
  
  prepared_data <- prepare_data(ps, group_var)
  raw_otu <- prepared_data$raw_otu
  sample_data_df <- prepared_data$sample_data_df
  y <- prepared_data$y
  
  raw_otu <- raw_otu + pseudocount
  lib.size <- colSums(raw_otu)
  scaling.factor <- lib.size / mean(lib.size)
  dat.normed <- t(t(raw_otu) / (scaling.factor + pseudocount))
  
  dge <- edgeR::DGEList(counts = raw_otu, group = y)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge)
  design <- stats::model.matrix(~ y)
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit)
  topTags <- edgeR::topTags(lrt, n = nrow(raw_otu))
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor, differential_abundance = topTags))
}

# -----------------------------------------------------------
#' TMM Normalization (Trimmed Mean of M component)
#' 
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @importFrom edgeR DGEList calcNormFactors
norm.TMM <- function(ps, groups) {
  otu_table_matrix <- as(phyloseq::otu_table(ps), "matrix")
  
  zero_rows <- rowSums(otu_table_matrix) == 0
  if (any(zero_rows)) {
    otu_table_matrix <- otu_table_matrix[!zero_rows, ]
  }
  
  zero_cols <- colSums(otu_table_matrix) == 0
  if (any(zero_cols)) {
    otu_table_matrix <- otu_table_matrix[, !zero_cols]
  }
  
  sample_data_df <- as(phyloseq::sample_data(ps), "data.frame")
  group <- sample_data_df[[groups]]
  dge <- edgeR::DGEList(counts = otu_table_matrix, group = group)
  
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  scaling.factor <- dge$samples$norm.factors
  dat.normed <- t(t(dge$counts) / scaling.factor)
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = phyloseq::taxa_are_rows(ps))
  phyloseq::sample_data(ps)$norm_factors <- scaling.factor
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' CLR Normalization (Centered Log-Ratio Transformation)
#' 
#' @param ps A phyloseq object.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.clr <- function(ps) {
  ps <- remove_zero_negative_count_samples(ps)
  
  gm_mean <- function(x) exp(sum(log(x[x > 0])) / length(x))
  
  ps_clr <- phyloseq::transform_sample_counts(ps, function(x) log(x / gm_mean(x)))
  scaling.factor <- rep(1, phyloseq::nsamples(ps_clr))
  
  return(list(dat.normed = ps_clr, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' Rarefying
#' 
#' @param ps A phyloseq object.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @importFrom phyloseq rarefy_even_depth sample_sums
norm.rar <- function(ps) {
  ps <- remove_zero_negative_count_samples(ps)
  
  ps_rarefied <- phyloseq::rarefy_even_depth(ps, rngseed = 123)
  scaling.factor <- phyloseq::sample_sums(ps_rarefied)
  
  return(list(dat.normed = ps_rarefied, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' TSS Normalization (Total Sum Scaling)
#' 
#' @param ps A phyloseq object.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.tss <- function(ps) {
  ps <- remove_zero_negative_count_samples(ps)
  
  otu <- phyloseq::otu_table(ps)
  size <- colSums(otu)
  otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
  phyloseq::otu_table(ps) <- phyloseq::otu_table(otu_normed, taxa_are_rows = phyloseq::taxa_are_rows(ps))
  
  scaling.factor <- rep(1, phyloseq::nsamples(ps))
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' CSS Normalization (Cumulative Sum Scaling)
#' 
#' @param ps A phyloseq object.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @importFrom edgeR DGEList calcNormFactors
norm.css <- function(ps) {
  ps <- remove_zero_negative_count_samples(ps)
  
  raw <- as(phyloseq::otu_table(ps), "matrix") + 1e-6
  dge <- edgeR::DGEList(counts = raw)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  scaling.factor <- dge$samples$norm.factors
  dat.normed <- t(t(raw) / scaling.factor)
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' RLE Normalization (Relative Log Expression)
#' 
#' @param ps A phyloseq object.
#' @param locfunc A function to compute the location statistic (default is median).
#' @param type A character string specifying the type of normalization ("poscounts" or "ratio").
#' @param geo_means A vector of geometric means for each feature.
#' @param control_genes A vector of control genes.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
norm.rle <- function(ps, locfunc = stats::median, type = c("poscounts", "ratio"), geo_means = NULL, control_genes = NULL) {
  type <- match.arg(type, c("poscounts", "ratio"))
  
  otu <- as(phyloseq::otu_table(ps), "matrix")
  
  if (is.null(geo_means)) {
    geo_means <- apply(otu, 1, gm_mean)
  }
  
  nf <- DESeq2::estimateSizeFactorsForMatrix(
    otu,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )
  
  phyloseq::otu_table(ps) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  ps <- set_nf(ps, nf)
  
  return(list(dat.normed = ps, scaling.factor = nf))
}

# -----------------------------------------------------------

# Example usage for TC normalization
# ps is phyloseq object and sample_data(ps)$Animal.type contains your group labels
# ps=physeq_16SOTU
# Host.species <- as.factor(ps@sam_data$Host.species)
# result_TC <- normalization_set(ps, method = "TC", groups = "Host.species")
# normalized_ps_TC <- result_TC$dat.normed
# scaling_factors_TC <- result_TC$scaling.factor

# -----------------------------------------------------------
# Example for UQ normalization
# result_UQ <- normalization_set(ps, method = "UQ", groups = "Host.species")
# normalized_ps_UQ <- result_UQ$dat.normed
# scaling_factors_UQ <- result_UQ$scaling.factor

# -----------------------------------------------------------
# Example for Median normalization
# result_med <- normalization_set(ps, method = "med", groups = "Host.species")
# normalized_ps_med <- result_med$dat.normed
# scaling_factors_med <- result_med$scaling.factor

# -----------------------------------------------------------
# Example for DESeq normalization
# ps_n<-remove_zero_negative_count_samples(ps)
# result_DESeq <- normalization_set(ps_n, method = "DESeq", groups = "Animal.type")
# normalized_ps_DESeq <- result_DESeq$dat.normed
# scaling_factors_DESeq <- result_DESeq$scaling.factor

# -----------------------------------------------------------
# Example for Poisson normalization
# result_Poisson <- normalization_set(ps, method = "Poisson", groups = "Host.genus")
# normalized_ps_Poisson <- result_Poisson$dat.normed
# scaling_factors_Poisson <- result_Poisson$scaling.factor

# -----------------------------------------------------------
# Example for Quantile normalization
# result_QN <- normalization_set(ps, method = "QN")
# normalized_ps_QN <- result_QN$dat.normed
# scaling_factors_QN <- result_QN$scaling.factor

# -----------------------------------------------------------
# Example for TMM normalization
# result_TMM <- normalization_set(physeq_ITS_adj_scaled_n, method = "TMM", groups = "Animal.type")
# normalized_ps_TMM <- result_TMM$dat.normed
# scaling_factors_TMM <- result_TMM$scaling.factor

# -----------------------------------------------------------
# Example for CLR normalization
# result_clr <- normalization_set(ps, method = "clr")
# normalized_ps_clr <- result_clr$dat.normed
# scaling_factors_clr <- result_clr$scaling.factor

# -----------------------------------------------------------
# Example for Rarefying
# result_rar <- normalization_set(ps, method = "rar")
# normalized_ps_rar <- result_rar$dat.normed
# scaling_factors_rar <- result_rar$scaling.factor

# -----------------------------------------------------------
# Example for CSS normalization
# result_css <- normalization_set(ps, method = "css")
# normalized_ps_css <- result_css$dat.normed
# scaling_factors_css <- result_css$scaling.factor

# -----------------------------------------------------------
# Example for TSS normalization
# result_tss <- normalization_set(ps, method = "tss")
# normalized_ps_tss <- result_tss$dat.normed
# scaling_factors_tss <- result_tss$scaling.factor

# -----------------------------------------------------------
# Example for RLE normalization
# result_rle <- normalization_set(ps, method = "rle")
# normalized_ps_rle <- result_rle$dat.normed
# scaling_factors_rle <- result_rle$scaling.factor

# -----------------------------------------------------------
