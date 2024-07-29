# For Differential Abundance Analysis: Consider DESeq, TMM, or CSS normalization, as these methods are designed to handle compositional biases and library size differences.
# For Compositional Data: CLR normalization is a good choice, as it transforms the data to account for its compositional nature.
# For Simplicity and Ease of Use: TC, UQ, or Median normalization are quick and straightforward but may not be as robust.
# For Handling Batch Effects: SVA or RUV methods are effective.

#' Install and Load Required Packages
#'
#' @param packages A character vector of package names to install and load.
install_and_load <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("phyloseq", "DESeq2", "edgeR", "sva", "EDASeq", "RUVSeq", "BiocManager", "BiocGenerics")) {
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

required_packages <- c("phyloseq", "DESeq2", "edgeR", "sva", "EDASeq", "RUVSeq", "BiocManager", "BiocGenerics")
install_and_load(required_packages)

#' Calculate Geometric Mean
#'
#' @param x A numeric vector.
#' @param na.rm Logical, should missing values (NAs) be removed?
#' @return Geometric mean of x.
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
}

#' Remove Samples with Zero, Negative Counts, or NA Values and Add Pseudocount
#'
#' @param ps A phyloseq object.
#' @param pseudocount A numeric value to add to avoid zero counts.
#' @return A phyloseq object with filtered and adjusted OTU table.
remove_zero_negative_count_samples <- function(ps, pseudocount = 1e-6) {
  otu <- as(otu_table(ps), "matrix")
  
  # Remove samples with zero, negative counts, or NA values
  zero_negative_count_samples <- sample_sums(ps) <= 0
  na_count_samples <- apply(otu, 2, function(x) any(is.na(x)))
  samples_to_remove <- zero_negative_count_samples | na_count_samples
  if (any(samples_to_remove)) {
    cat("Removing", sum(samples_to_remove), "samples with zero, negative counts, or NA values.\n")
    ps <- prune_samples(!samples_to_remove, ps)
    otu <- as(otu_table(ps), "matrix")
  }
  
  # Remove rows (features) that are completely zero across all samples
  zero_rows <- rowSums(otu) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples.\n")
    otu <- otu[!zero_rows, ]
    ps <- prune_taxa(!zero_rows, ps)
  }
  
  # Add pseudocount to avoid zero counts
  otu <- otu + pseudocount
  
  # Round counts to ensure they are integers
  otu <- round(otu)
  
  otu_table(ps) <- otu_table(otu, taxa_are_rows = TRUE)
  
  return(ps)
}

#' Convert Categorical Columns to Factors in Sample Data
#'
#' @param ps A phyloseq object.
#' @return A phyloseq object with updated sample data.
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

#' Create a List from a Phyloseq Object
#'
#' @param physeq A phyloseq object.
#' @return A list containing the DGE list and updated phyloseq object.
create_list <- function(physeq) {
  if (!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }
  counts <- as(otu_table(physeq), "matrix")
  sample_data_df <- as(sample_data(physeq), "data.frame")
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
  sample_data(updated_physeq) <- sample_data(sample_data_df)
  return(list(dge_list = dge_base, phyloseq_obj = updated_physeq))
}

#' Apply the Selected Normalization Method to the Phyloseq Object
#'
#' @param ps A phyloseq object.
#' @param method A character string specifying the normalization method ("TC", "UQ", "med", "DESeq", "Poisson", "QN", "SVA", "RUVg", "RUVs", "RUVr", "TMM", "clr", "rar", "css", "tss", "rle").
#' @param groups A column name of group labels from sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @examples
#' result_TMM <- normalization_set(physeq_16SASV, method = "TMM", groups = "Animal.type")
#' normalized_ps_TMM <- result_TMM$dat.normed
#' scaling_factors_TMM <- result_TMM$scaling.factor
#' @export
normalization_set <- function(ps, method, groups = NULL) {
  if (nsamples(ps) == 0) stop("The phyloseq object contains no samples.")
  if (length(sample_names(ps)) == 0) stop("Sample names are missing in the phyloseq object.")
  
  ps <- remove_zero_negative_count_samples(ps)
  ps <- convert_categorical_to_factors(ps)
  
  # Ensure groups is a single column name from sample data
  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) != 1 || !groups %in% colnames(sample_data(ps))) {
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
                   "SVA" = norm.SVA(ps, groups),
                   "RUVg" = norm.RUVg(ps, groups),
                   "RUVs" = norm.RUVs(ps, groups),
                   "RUVr" = norm.RUVr(ps, groups),
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
  
  # Ensure the scaling factors are named and match the sample names in the phyloseq object
  if (length(scaling.factor) != nsamples(ps)) {
    stop("Length of scaling.factor must match the number of samples in the phyloseq object.")
  }
  
  # Add the scaling factors to the sample data
  sample_data(ps)$norm_factors <- scaling.factor
  return(ps)
}

#' Process Data with Feature Category
#'
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all_features", "interquartile_range", "non_zero_counts", "specified_features").
#' @return A processed phyloseq object.
process_data_with_feature_category <- function(ps, feature_category = c("all_features", "interquartile_range", "non_zero_counts", "specified_features")) {
  feature_category <- match.arg(feature_category)
  if (feature_category == "all_features") {
    return(ps)
  }
  # Additional filtering criteria can be added here based on the selected feature category.
  stop("Feature category processing not implemented.")
}

#' Poisson Normalization and Differential Abundance Function
#'
#' @param ps A phyloseq object or matrix of raw counts.
#' @param group_var A string specifying the grouping variable in sample data (if phyloseq object).
#' @param pseudocount A numeric value added to avoid division by zero.
#' @return A list containing the normalized data, scaling factor, and differential abundance results.
#' @export
norm.Poisson <- function(ps, group_var = NULL, pseudocount = 1e-6) {
  ps <- remove_zero_negative_count_samples(ps)
  
  prepare_data <- function(ps, group_var) {
    if (inherits(ps, "phyloseq")) {
      raw_otu <- as(otu_table(ps), "matrix")
      sample_data_df <- as(sample_data(ps), "data.frame")
      
      if (!is.null(group_var)) {
        if (!is.character(group_var) || length(group_var) != 1 || !group_var %in% colnames(sample_data_df)) {
          stop("The specified group_var does not exist in the sample data.")
        }
        y <- as.numeric(as.factor(sample_data_df[[group_var]]))
      } else {
        y <- rep(1, ncol(raw_otu))  # All samples in one group if group_var is NULL
      }
    } else if (is.matrix(ps)) {
      raw_otu <- ps
      if (!is.null(group_var)) {
        stop("group_var should be NULL when ps is a matrix.")
      }
      y <- rep(1, ncol(raw_otu))  # All samples in one group
    } else {
      stop("Input ps must be a phyloseq object or a matrix.")
    }
    list(raw_otu = raw_otu, sample_data_df = sample_data_df, y = y)
  }
  
  prepared_data <- prepare_data(ps, group_var)
  raw_otu <- prepared_data$raw_otu
  sample_data_df <- prepared_data$sample_data_df
  y <- prepared_data$y
  
  # Add pseudocount again to avoid division by zero
  raw_otu <- raw_otu + pseudocount
  
  lib.size <- colSums(raw_otu)
  scaling.factor <- lib.size / mean(lib.size)
  dat.normed <- t(t(raw_otu) / (scaling.factor + pseudocount))
  dge_base <- list(
    counts = raw_otu,
    samples = sample_data_df,
    genes = NULL,
    group = if (!is.null(group_var)) sample_data_df[[group_var]] else NULL,
    lib.size = lib.size,
    norm.factors = scaling.factor
  )
  dge <- DGEList(counts = dge_base$counts, group = dge_base$group)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge)
  design <- model.matrix(~ y)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit)
  topTags <- topTags(lrt, n = nrow(raw_otu))
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor, differential_abundance = topTags))
}

#' TMM Normalization (Trimmed Mean of M component)
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references Robinson and Oshlack, 2010.
#' @export
norm.TMM <- function(ps, groups) {
  # Ensure the OTU table is numeric
  otu_table_matrix <- as(otu_table(ps), "matrix")
  if (!is.numeric(otu_table_matrix)) {
    otu_table_matrix <- apply(otu_table_matrix, 2, as.numeric)
    if (any(is.na(otu_table_matrix))) {
      stop("OTU table contains non-numeric values that could not be converted.")
    }
  }
  
  # Remove rows (features) that are completely zero across all samples
  zero_rows <- rowSums(otu_table_matrix) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples.\n")
    otu_table_matrix <- otu_table_matrix[!zero_rows, ]
  }
  
  # Remove columns (samples) that are completely zero across all features
  zero_cols <- colSums(otu_table_matrix) == 0
  if (any(zero_cols)) {
    cat("Removing", sum(zero_cols), "samples with zero counts across all features.\n")
    otu_table_matrix <- otu_table_matrix[, !zero_cols]
  }
  
  # Check for extreme values
  if (any(is.infinite(otu_table_matrix)) || any(is.nan(otu_table_matrix))) {
    stop("OTU table contains infinite or NaN values. Please clean your data.")
  }
  
  # Extract sample data
  sample_data_df <- as(sample_data(ps), "data.frame")
  
  # Ensure the group variable is specified correctly
  if (!groups %in% colnames(sample_data_df)) {
    stop("The specified group variable does not exist in the sample data.")
  }
  
  # Align the sample data with the filtered OTU table
  sample_data_df <- sample_data_df[colnames(otu_table_matrix), , drop = FALSE]
  
  # Prepare the DGEList object
  group <- sample_data_df[[groups]]
  dge <- DGEList(counts = otu_table_matrix, group = group)
  
  # Perform TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Extract the normalization factors
  scaling.factor <- dge$samples$norm.factors
  
  # Normalize the OTU table
  dat.normed <- t(t(dge$counts) / scaling.factor)
  
  # Determine if taxa are rows
  taxa_are_rows_flag <- taxa_are_rows(otu_table(ps))
  
  # Create a new OTU table preserving the taxa_are_rows attribute
  otu_table_norm <- otu_table(dat.normed, taxa_are_rows = taxa_are_rows_flag)
  
  # Update the phyloseq object with normalized data
  otu_table(ps) <- otu_table_norm
  sample_data(ps)$norm.factors <- scaling.factor
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' TC Normalization (Total Count Scaling)
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @export
norm.TC <- function(ps, groups) {
  ps <- create_list(ps)$phyloseq_obj
  dat.DGE <- create_list(ps)$dge_list
  scaling.factor <- dat.DGE$samples$lib.size / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' UQ Normalization (Upper Quartile)
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @export
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
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  
  # Set normalization factors
  ps <- set_nf(ps, scaling.factor)
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' Median Normalization
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @export
norm.med <- function(ps, groups) {
  ps <- create_list(ps)$phyloseq_obj
  dat.DGE <- create_list(ps)$dge_list
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' DESeq Normalization
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references DESeq2 package.
#' @export
norm.DESeq <- function(ps, groups) {
  ps <- remove_zero_negative_count_samples(ps, pseudocount = 1e-6)  # Ensure no zero counts and add pseudocount
  raw <- as(otu_table(ps), "matrix")
  
  # Check if the groups vector has only one unique value
  unique_groups <- unique(groups)
  if (length(unique_groups) == 1) {
    # All samples have the same condition, use design ~ 1
    design <- ~ 1
  } else {
    design <- ~ Condition
  }
  
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) <- colnames(raw)
  dat.DGE <- DESeqDataSetFromMatrix(countData = raw, colData = condition, design = design)
  
  # Handle case when all samples have the same condition
  if (length(unique_groups) == 1) {
    dat.DGE <- estimateSizeFactors(dat.DGE)
  } else {
    dat.DGE <- DESeq(dat.DGE, fitType = "local")
  }
  
  scaling.factor <- sizeFactors(dat.DGE)
  dat.normed <- counts(dat.DGE, normalized = TRUE)
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' Quantile Normalization (QN) for phyloseq object using base R
#'
#' @param ps A phyloseq object.
#' @param filter Logical, whether to filter low counts.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @export
norm.QN <- function(ps, filter = FALSE) {
  otu <- as(otu_table(ps), "matrix")
  
  if (!is.numeric(otu)) {
    stop("OTU table must contain numeric values.")
  }
  
  taxa_names_original <- taxa_names(ps)
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
  otu_table(ps) <- otu_table(normalized, taxa_are_rows = TRUE)
  
  # for QN we don't need scaling factors
  scaling.factor <- NULL
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' Surrogate Variable Analysis for Sequencing Data (SVA)
#'
#' @param ps A phyloseq object.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and adjustment factors.
#' @references SVA package.
#' @export
norm.SVA <- function(ps, groups) {
  # Ensure 'groups' is a valid character vector
  if (is.null(groups) || !is.character(groups)) {
    stop("'groups' must be a character vector specifying a single column name from sample data.")
  }
  
  # Convert groups to a factor
  group_var <- sample_data(ps)[[groups]]
  if (!is.factor(group_var)) {
    group_var <- as.factor(group_var)
  }
  
  # Debugging information
  cat("Levels in 'groups':", levels(group_var), "\n")
  cat("Length of 'group_var':", length(group_var), "\n")
  
  # Check the number of levels in the groups factor
  if (length(levels(group_var)) < 2) {
    stop("The 'groups' factor must have at least two levels.")
  }
  
  raw <- as(otu_table(ps), "matrix")
  
  # Ensure the number of columns in raw matches the number of samples
  if (ncol(raw) != nsamples(ps)) {
    stop("The number of columns in 'raw' does not match the number of samples in the phyloseq object.")
  }
  
  # Add a small constant to all counts to avoid zero quantiles
  raw <- raw + 1e-6
  
  # Additional filtering criteria
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.sva <- raw[filter, ]
  
  # Align the sample data with the filtered OTU table
  valid_samples <- colnames(dat.sva)
  group_var <- group_var[valid_samples]
  sample_data_filtered <- sample_data(ps)[valid_samples, , drop = FALSE]
  
  # Ensure the number of columns in dat.sva matches the number of samples in group_var
  dat.sva <- dat.sva[, rownames(sample_data_filtered)]
  
  # Debugging information
  cat("Dimensions of filtered 'dat.sva':", dim(dat.sva), "\n")
  cat("Length of filtered 'group_var':", length(group_var), "\n")
  
  mod1 <- model.matrix(~ group_var)
  mod0 <- cbind(mod1[,1])
  dat0 <- as.matrix(dat.sva)
  
  # Ensure the dimensions match
  if (ncol(dat0) != nrow(mod1)) {
    stop("The number of columns in 'dat0' does not match the number of rows in 'mod1'.")
  }
  
  # Catch errors when running svaseq and provide a meaningful error message
  svseq <- tryCatch(
    {
      invisible(capture.output(svaseq(dat0, mod1, mod0, n.sv = 1)$sv))
    },
    error = function(e) {
      stop("SVA analysis failed: ", e$message)
    }
  )
  
  # More debugging information
  cat("Dimensions of 'svseq':", dim(svseq), "\n")
  
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(raw))
  P <- ncol(mod1)
  dat.normed <- raw - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  
  # Create the phyloseq object
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  
  # Remove rows with zero counts from the phyloseq object
  zero_rows <- rowSums(otu_table(ps)) == 0
  if (any(zero_rows)) {
    cat("Removing", sum(zero_rows), "features with zero counts across all samples after SVA transformation.\n")
    ps <- prune_taxa(!zero_rows, ps)
  }
  
  # Assuming set_nf is a placeholder for any further normalization steps or attributes setting
  ps <- set_nf(ps, svseq)
  
  return(list(dat.normed = ps, adjust.factor = svseq))
}

#' Remove Unwanted Variation Using Control Genes (RUVg)
#'
#' @param ps A phyloseq object.
#' @param groups The name of the sample variable in the sample data that indicates group labels.
#' @return A list containing the normalized phyloseq object and adjustment factors.
#' @references Gagnon-Bartsch, Jacob, and Speed 2013; Risso et al. 2014; Gagnon-Bartsch and Speed 2012.
#' @export
norm.RUVg <- function(ps, groups) {
  # Remove samples with zero, negative counts, or NA values
  ps <- remove_zero_negative_count_samples(ps)
  
  # Extract OTU table
  raw <- as(otu_table(ps), "matrix")
  print("Dimensions of raw OTU table:")
  print(dim(raw))
  
  # Apply filter to remove features with zero counts across all samples
  raw <- raw[rowSums(raw) > 0, ]
  print("Dimensions of OTU table after removing features with zero counts:")
  print(dim(raw))
  
  # Apply filter to remove samples with zero counts across all features
  raw <- raw[, colSums(raw) > 0]
  print("Dimensions of OTU table after removing samples with zero counts:")
  print(dim(raw))
  
  # Add a small constant to all counts to avoid zero quantiles
  raw <- raw + 1e-6
  
  # Apply additional filtering criteria
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  print("Dimensions of filtered OTU table:")
  print(dim(dat.ruv))
  
  # Extract groups from sample data
  condition <- factor(sample_data(ps)[[groups]])
  print("Length of condition vector:")
  print(length(condition))
  
  # Check if the number of columns in dat.ruv matches the length of condition
  if (ncol(dat.ruv) != length(condition)) {
    stop("The length of 'groups' must match the number of columns in the OTU table.")
  }
  
  # Create SeqExpressionSet
  set <- newSeqExpressionSet(as.matrix(dat.ruv), phenoData = data.frame(condition, row.names = colnames(dat.ruv)))
  print("Dimensions of SeqExpressionSet:")
  print(dim(counts(set)))
  
  design <- model.matrix(~ condition, data = data.frame(condition, row.names = colnames(dat.ruv)))
  print("Design matrix:")
  print(design)
  
  y <- DGEList(counts = as.matrix(counts(set)), group = condition)
  
  # Remove samples with zero or NA library sizes
  zero_lib_size <- y$samples$lib.size == 0 | is.na(y$samples$lib.size)
  if (any(zero_lib_size)) {
    cat("Removing", sum(zero_lib_size), "samples with zero or NA library sizes.\n")
    y <- y[, !zero_lib_size]
  }
  
  # Check if all library sizes are finite
  if (any(!is.finite(y$samples$lib.size))) {
    stop("Library sizes must be finite.")
  }
  
  # Use the TMM method for normalization
  y <- calcNormFactors(y, method = "TMM")
  
  # Ensure all normalization factors are finite
  if (any(!is.finite(y$samples$norm.factors))) {
    stop("Normalization factors must be finite.")
  }
  
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1/(0.15 * nrow(raw))]))]
  t <- RUVg(set, spikes, k = 1)
  dat.normed <- normCounts(t)
  
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, t$W)
  
  return(list(dat.normed = ps, adjust.factor = t$W))
}

#' Remove Unwanted Variation Using Replicate Samples (RUVs)
#'
#' @param ps A phyloseq object.
#' @param groups A vector of group labels.
#' @return A list containing the normalized phyloseq object and adjustment factors.
#' @references Gagnon-Bartsch, Jacob, and Speed 2013; Risso et al. 2014; Gagnon-Bartsch and Speed 2012.
#' @export
norm.RUVs <- function(ps, groups) {
  print("Starting norm.RUVr function")
  
  raw <- as(otu_table(ps), "matrix")
  print("OTU table dimensions:")
  print(dim(raw))
  
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  print("Filtered OTU table dimensions:")
  print(dim(dat.ruv))
  
  genes <- rownames(dat.ruv)
  condition <- sample_data(ps)[[groups]]
  print("Condition vector length:")
  print(length(condition))
  
  set <- newSeqExpressionSet(as.matrix(dat.ruv), phenoData = AnnotatedDataFrame(data.frame(condition, row.names = colnames(dat.ruv))))
  print("SeqExpressionSet dimensions:")
  print(dim(counts(set)))
  
  design <- model.matrix(~ condition, data = pData(set))
  print("Design matrix:")
  print(design)
  
  y <- DGEList(counts = as.matrix(counts(set)), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type = "deviance")
  setUQ <- betweenLaneNormalization(set, which = "upper")
  controls <- rownames(dat.ruv)
  t <- RUVr(setUQ, controls, k = 1, res)
  dat.normed <- normCounts(t)
  
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, t$W)
  
  print("Completed norm.RUVr function")
  return(list(dat.normed = ps, adjust.factor = t$W))
}

#' CLR Normalization (Centered Log-Ratio Transformation)
#'
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all_features", "interquartile_range", "non_zero_counts", "specified_features").
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references Aitchison, 1986.
#' @export
norm.clr <- function(ps, feature_category = c("all_features", "interquartile_range", "non_zero_counts", "specified_features")) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category)
  if (nrow(otu_table(ps)) == 0 || ncol(otu_table(ps)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  gm_mean <- function(x) exp(sum(log(x[x > 0])) / length(x))
  ps_clr <- transform_sample_counts(ps, function(x) log(x / gm_mean(x)))
  scaling.factor <- rep(1, nsamples(ps_clr))
  return(list(dat.normed = ps_clr, scaling.factor = scaling.factor))
}

#' Rarefying
#'
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all_features", "interquartile_range", "non_zero_counts", "specified_features")).
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references McMurdie and Holmes, 2014.
#' @export
norm.rar <- function(ps, feature_category = c("all_features", "interquartile_range", "non_zero_counts", "specified_features")) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category)
  if (nrow(otu_table(ps)) == 0 || ncol(otu_table(ps)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  cat("Performing rarefaction with set.seed(123) for reproducibility.\n")
  ps_rarefied <- rarefy_even_depth(ps, rngseed = 123)
  scaling.factor <- sample_sums(ps_rarefied)
  return(list(dat.normed = ps_rarefied, scaling.factor = scaling.factor))
}

#' TSS Normalization (Total Sum Scaling)
#'
#' @param ps A phyloseq object.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references Bolstad et al., 2003.
#' @export
norm.tss <- function(ps) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category = "all_features")
  
  otu <- otu_table(ps)
  size <- colSums(otu)
  otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
  otu_table(ps) <- otu_table(otu_normed, taxa_are_rows = taxa_are_rows(ps))
  scaling.factor <- rep(1, nsamples(ps))  # Since TSS normalizes to relative abundance
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' CSS Normalization (Cumulative Sum Scaling)
#'
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all_features", "interquartile_range", "non_zero_counts", "specified_features").
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references Paulson et al., 2013.
#' @export
norm.css <- function(ps, feature_category = c("all_features", "interquartile_range", "non_zero_counts", "specified_features")) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category)
  
  if (nrow(otu_table(ps)) == 0 || ncol(otu_table(ps)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  
  cat("Initial OTU table dimensions: ", dim(otu_table(ps)), "\n")
  raw <- as(otu_table(ps), "matrix")
  raw <- raw + 1e-6  # Add pseudocount
  dat.DGE <- DGEList(counts = raw)
  dat.DGE <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- dat.DGE$samples$norm.factors
  dat.normed <- t(t(raw) / scaling.factor)
  
  # Ensure dat.normed is converted back to otu_table format
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  ps <- set_nf(ps, scaling.factor)
  
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' RLE Normalization (Relative Log Expression)
#'
#' @param ps A phyloseq object.
#' @param locfunc A function to compute the location statistic (default is median).
#' @param type A character string specifying the type of normalization ("poscounts" or "ratio").
#' @param geo_means A vector of geometric means for each feature.
#' @param control_genes A vector of control genes.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @references Anders and Huber, 2010.
#' @export
norm.rle <- function(ps, locfunc = stats::median, type = c("poscounts", "ratio"), geo_means = NULL, control_genes = NULL) {
  type <- match.arg(type, c("poscounts", "ratio"))
  otu <- as(otu_table(ps), "matrix")
  
  if (is.null(geo_means)) {
    geo_means <- apply(otu, 1, gm_mean)
  }
  
  # Ensure control_genes is defined
  if (is.null(control_genes)) {
    control_genes <- rep(TRUE, nrow(otu))  # Use all genes as control by default
  } else if (is.numeric(control_genes)) {
    control_genes <- control_genes <= nrow(otu)
  } else if (is.logical(control_genes)) {
    if (length(control_genes) != nrow(otu)) {
      stop("control_genes should be either a numeric vector of indices or a logical vector of length equal to the number of rows of counts.")
    }
  } else {
    stop("control_genes should be either a numeric or logical vector.")
  }
  
  nf <- estimateSizeFactorsForMatrix(
    otu,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )
  
  otu_table(ps) <- otu_table(otu, taxa_are_rows = TRUE)
  ps <- set_nf(ps, nf)
  return(list(dat.normed = ps, scaling.factor = nf))
}


# # Example usage for TC normalization
# # Assuming ps is your phyloseq object and sample_data(ps)$Animal.type contains your group labels
# ps <- SP_Plethodon
# Host.species <- as.factor(ps@sam_data$Host.species)
# ps@sam_data
# result_TC <- normalization_set(ps, method = "TC", groups = "Host.species")
# normalized_ps_TC <- result_TC$dat.normed
# scaling_factors_TC <- result_TC$scaling.factor
# 
# # Example for UQ normalization
# result_UQ <- normalization_set(ps, method = "UQ", groups = "Host.species")
# normalized_ps_UQ <- result_UQ$dat.normed
# scaling_factors_UQ <- result_UQ$scaling.factor
# 
# # Example for Median normalization
# result_med <- normalization_set(ps, method = "med", groups = "Host.species")
# normalized_ps_med <- result_med$dat.normed
# scaling_factors_med <- result_med$scaling.factor
# 
# # Example for DESeq normalization
# result_DESeq <- normalization_set(ps, method = "DESeq", groups = "Host.species")
# normalized_ps_DESeq <- result_DESeq$dat.normed
# scaling_factors_DESeq <- result_DESeq$scaling.factor
# 
# # Example for Poisson normalization
# result_Poisson <- normalization_set(ps, method = "Poisson", groups = "Host.species")
# normalized_ps_Poisson <- result_Poisson$dat.normed
# scaling_factors_Poisson <- result_Poisson$scaling.factor
# 
# # Example for Quantile normalization
# result_QN <- normalization_set(ps, method = "QN")
# normalized_ps_QN <- result_QN$dat.normed
# scaling_factors_QN <- result_QN$scaling.factor
# 
# # Example for SVA normalization
# result_SVA <- normalization_set(ps, method = "SVA", groups = "Host.species")
# normalized_ps_SVA <- result_SVA$dat.normed
# scaling_factors_SVA <- result_SVA$scaling.factor
# 
# # Example for RUVg normalization
# result_RUVg <- normalization_set(ps, method = "RUVg", groups = "Host.species")
# normalized_ps_RUVg <- result_RUVg$dat.normed
# scaling_factors_RUVg <- result_RUVg$scaling.factor
# 
# # Example for RUVs normalization
# result_RUVs <- normalization_set(ps, method = "RUVs", groups = "Host.species")
# normalized_ps_RUVs <- result_RUVs$dat.normed
# scaling_factors_RUVs <- result_RUVs$scaling.factor
# 
# # Example for RUVr normalization
# result_RUVr <- normalization_set(ps, method = "RUVr", groups = "Host.species")
# normalized_ps_RUVr <- result_RUVr$dat.normed
# scaling_factors_RUVr <- result_RUVr$scaling.factor
# 
# # Example for TMM normalization
# result_TMM <- normalization_set(physeq_16SASV, method = "TMM", groups = "Animal.type")
# normalized_ps_TMM <- result_TMM$dat.normed
# scaling_factors_TMM <- result_TMM$scaling.factor
# 
# # Example for CLR normalization
# result_clr <- normalization_set(ps, method = "clr")
# normalized_ps_clr <- result_clr$dat.normed
# scaling_factors_clr <- result_clr$scaling.factor
# 
# # Example for Rarefying
# result_rar <- normalization_set(ps, method = "rar")
# normalized_ps_rar <- result_rar$dat.normed
# scaling_factors_rar <- result_rar$scaling.factor
# 
# # Example for CSS normalization
# result_css <- normalization_set(ps, method = "css")
# normalized_ps_css <- result_css$dat.normed
# scaling_factors_css <- result_css$scaling.factor
# 
# # Example for TSS normalization
# result_tss <- normalization_set(ps, method = "tss")
# normalized_ps_tss <- result_tss$dat.normed
# scaling_factors_tss <- result_tss$scaling.factor
# 
# # Example for RLE normalization
# result_rle <- normalization_set(ps, method = "rle")
# normalized_ps_rle <- result_rle$dat.normed
# scaling_factors_rle <- result_rle$scaling.factor
# 
