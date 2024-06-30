# Function to install and load required packages
install_and_load <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (package %in% c("phyloseq", "DESeq2", "edgeR", "PoissonSeq", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics")) {
        BiocManager::install(package)
      } else {
        install.packages(package)
      }
    }
    library(package, character.only = TRUE)
  }
}

# Specific installation for PoissonSeq
if (!requireNamespace("PoissonSeq", quietly = TRUE)) {
  install.packages("https://cran.r-project.org/src/contrib/PoissonSeq_1.1.2.tar.gz", repos = NULL, type = "source")
}

# List of required packages
required_packages <- c("phyloseq", "DESeq2", "edgeR", "PoissonSeq", "preprocessCore", "sva", "EDASeq", "RUVSeq", "Biobase", "BiocGenerics", "vegan", "chemometrics")

# Install and load required packages
install_and_load(required_packages)

# Helper function for geometric mean
#' Helper function for geometric mean
#' @param x A numeric vector.
#' @param na.rm Logical, whether to remove NA values.
#' @return Geometric mean of the input vector.
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm = na.rm) / length(x))
}

# Remove samples with zero or negative counts
#' Remove samples with zero or negative counts
#' @param ps A phyloseq object.
#' @return A pruned phyloseq object without samples with zero or negative counts.
remove_zero_negative_count_samples <- function(ps) {
  zero_negative_count_samples <- sample_sums(ps) <= 0
  if (any(zero_negative_count_samples)) {
    cat("Removing", sum(zero_negative_count_samples), "samples with zero or negative counts.\n")
    ps <- prune_samples(!zero_negative_count_samples, ps)
  }
  return(ps)
}

# Process data with various feature categories
#' Process data with various feature categories
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all", "iqlr", "zero", "lvha").
#' @return A pruned phyloseq object based on the specified feature category.
process_data_with_feature_category <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha")) {
  feature_category <- match.arg(feature_category)
  
  if (feature_category == "all") {
    return(ps)
  } else if (feature_category == "iqlr") {
    variances <- apply(otu_table(ps), 1, var)
    lower_quartile <- quantile(variances, 0.25)
    upper_quartile <- quantile(variances, 0.75)
    keep_features <- variances >= lower_quartile & variances <= upper_quartile
    return(prune_taxa(keep_features, ps))
  } else if (feature_category == "zero") {
    counts <- otu_table(ps)
    keep_features <- apply(counts, 1, function(x) sum(x > 0) > 0)
    return(prune_taxa(keep_features, ps))
  } else if (feature_category == "lvha") {
    specified_taxa <- c("taxa1", "taxa2", "taxa3")
    keep_features <- rownames(otu_table(ps)) %in% specified_taxa
    return(prune_taxa(keep_features, ps))
  } else {
    stop("Invalid 'feature_category' option specified.")
  }
}

# Normalization set function
#' Normalization set function
#' @description Apply the selected normalization method to the phyloseq object.
#' @param ps A phyloseq object.
#' @param method A character string specifying the normalization method ("TC", "UQ", "med", "DESeq", "PoissonSeq", "QN", "SVA", "RUVg", "RUVs", "RUVr", "TMM", "clr", "rar", "css", "tss", "rle").
#' @param groups A vector of group labels.
#' @return A list containing the normalized phyloseq object and scaling factors.
normalization_set <- function(ps, method, groups = NULL) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category = "all")
  raw <- as(otu_table(ps), "matrix")
  raw[is.na(raw)] <- 0  # Replace NA values with zero
  raw <- raw + 1  # Add pseudocount to all counts
  
  result <- switch(method,
                   "TC" = norm.TC(raw, groups),
                   "UQ" = norm.UQ(raw, groups),
                   "med" = norm.med(raw, groups),
                   "DESeq" = norm.DESeq(raw, groups),
                   "PoissonSeq" = norm.PoissonSeq(raw),
                   "QN" = norm.QN(raw),
                   "SVA" = norm.SVA(raw, groups),
                   "RUVg" = norm.RUVg(raw, groups),
                   "RUVs" = norm.RUVs(raw, groups),
                   "RUVr" = norm.RUVr(raw, groups),
                   "TMM" = norm.TMM(raw, groups),
                   "clr" = norm.clr(ps),
                   "rar" = norm.rar(ps),
                   "css" = norm.css(ps),
                   "tss" = norm.tss(ps),
                   "rle" = norm.rle(ps),
                   stop("Invalid normalization method"))
  
  dat.normed <- result$dat.normed
  scaling.factor <- result$scaling.factor
  
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  if (!is.null(scaling.factor)) {
    ps <- set_nf(ps, scaling.factor)
  }
  return(list(dat.normed = ps, scaling.factor = scaling.factor))
}

#' TMM Normalization (Trimmed Mean of M component)
#' @description A scaling normalization method used for RNA-seq count data to account for compositional differences between libraries.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and scaling factor.
#' @references Robinson and Oshlack, 2010.
norm.TMM <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)), group = factor(groups), genes = rownames(raw))
  dat.DGE <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- dat.DGE$samples$norm.factors
  dat.normed <- t(t(raw) / scaling.factor)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' TC Normalization (Total Count Scaling)
#' @description A simple scaling normalization method dividing by the library size.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and scaling factor.
norm.TC <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)), group = factor(groups), genes = rownames(raw))
  scaling.factor <- dat.DGE$samples$lib.size / 1e6
  dat.normed <- t(t(raw) / scaling.factor)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' UQ Normalization (Upper Quartile)
#' @description Normalizes data based on the upper quartile of the counts.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and scaling factor.
norm.UQ <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)), group = factor(groups), genes = rownames(raw))
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor / 1e6
  dat.normed <- t(t(raw) / scaling.factor)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' Median Normalization
#' @description Normalizes data based on the median count.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and scaling factor.
norm.med <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)), group = factor(groups), genes = rownames(raw))
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6
  dat.normed <- t(t(raw) / scaling.factor)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' DESeq Normalization
#' @description Normalizes counts based on the assumption that most genes are not differentially expressed.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and scaling factor.
#' @references DESeq2 package.
norm.DESeq <- function(raw, groups) {
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) <- colnames(raw)
  dat.DGE <- DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- estimateSizeFactors(dat.DGE)
  scaling.factor <- sizeFactors(dat.DGE)
  dat.normed <- counts(dat.DGE, normalized = TRUE)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' PoissonSeq Normalization
#' @description A normalization method specifically designed for RNA-seq data, available in the PoissonSeq package.
#' @param raw A matrix of raw counts.
#' @return A list containing the normalized data and scaling factor.
#' @references PoissonSeq Package.
norm.PoissonSeq <- function(raw) {
  invisible(capture.output(scaling.factor <- PS.Est.Depth(raw)))
  dat.normed <- t(t(raw) / scaling.factor)
  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}

#' Quantile Normalization (QN)
#' @description Adjusts the distribution of counts so that different samples have the same distribution.
#' @param raw A matrix of raw counts.
#' @param filter Logical, whether to filter low counts.
#' @return A list containing the normalized data.
#' @references Bolstad et al., 2003.
norm.QN <- function(raw, filter = FALSE) {
  if (filter == TRUE) {
    raw <- log2(raw + 1)
    raw <- raw[rowMeans(raw) > 2, ]
  } else {
    raw <- log2(raw + 1)
  }
  dat.log.normed <- normalize.quantiles(as.matrix(raw))
  dat.normed <- 2^dat.log.normed - 1
  colnames(dat.normed) <- colnames(raw)
  rownames(dat.normed) <- rownames(raw)
  return(list(dat.normed = dat.normed, scaling.factor = NULL))
}

#' Surrogate Variable Analysis for Sequencing Data (SVA)
#' @description Identifies and removes unwanted variation in sequencing data.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and adjustment factor.
#' @references SVA package.
norm.SVA <- function(raw, groups) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.sva <- raw[filter, ]
  genes <- rownames(dat.sva)
  mod1 <- model.matrix(~ groups)
  mod0 <- cbind(mod1[,1])
  dat0 <- as.matrix(dat.sva)
  invisible(capture.output(svseq <- svaseq(dat0, mod1, mod0, n.sv = 1)$sv))
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(raw))
  P <- ncol(mod1)
  dat.normed <- raw - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  return(list(dat.normed = dat.normed, adjust.factor = svseq))
}

#' Remove Unwanted Variation Using Control Genes (RUVg)
#' @description Uses control genes to remove unwanted variation.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and adjustment factor.
#' @references Gagnon-Bartsch, Jacob, and Speed 2013; Risso et al. 2014; Gagnon-Bartsch and Speed 2012.
norm.RUVg <- function(raw, groups) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv), phenoData = data.frame(condition, row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition, row.names = colnames(dat.ruv)))
  y <- DGEList(counts = as.matrix(counts(set)), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1/(0.15*nrow(raw))]))]
  t <- RUVg(set, spikes, k = 1)
  dat.normed <- normCounts(t)
  return(list(dat.normed = dat.normed, adjust.factor = t$W))
}

#' Remove Unwanted Variation Using Replicate Samples (RUVs)
#' @description Uses replicate samples to remove unwanted variation.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and adjustment factor.
#' @references Gagnon-Bartsch, Jacob, and Speed 2013; Risso et al. 2014; Gagnon-Bartsch and Speed 2012.
norm.RUVs <- function(raw, groups) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv), phenoData = data.frame(condition, row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition, row.names = colnames(dat.ruv)))
  y <- DGEList(counts = as.matrix(counts(set)), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1/(0.15*nrow(raw))]))]
  differences <- makeGroups(condition)
  controls <- rownames(dat.ruv)
  t <- RUVs(set, controls, k = 1, differences)
  dat.normed <- normCounts(t)
  return(list(dat.normed = dat.normed, adjust.factor = t$W))
}

#' Remove Unwanted Variation Using Residuals (RUVr)
#' @description Uses residuals to remove unwanted variation.
#' @param raw A matrix of raw counts.
#' @param groups A vector of group labels.
#' @return A list containing the normalized data and adjustment factor.
#' @references Gagnon-Bartsch, Jacob, and Speed 2013; Risso et al. 2014; Gagnon-Bartsch and Speed 2012.
norm.RUVr <- function(raw, groups) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv), phenoData = data.frame(condition, row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = pData(set))
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
  return(list(dat.normed = dat.normed, adjust.factor = t$W))
}

#' CLR Normalization (Centered Log-Ratio Transformation)
#' @description Computes log-ratios relative to the geometric mean of all features.
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all", "iqlr", "zero", "lvha").
#' @return A list containing the normalized data and scaling factor.
#' @references Aitchison, 1986.
norm.clr <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha")) {
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
#' @description Randomly removes reads from different samples until they all have the same predefined number of reads, ensuring equal library size.
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all", "iqlr", "zero", "lvha").
#' @return A list containing the normalized data and scaling factor.
#' @references McMurdie and Holmes, 2014.
norm.rar <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha")) {
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
#' @description Converts the feature table into relative abundance by dividing the total reads of each sample.
#' @param object A phyloseq object.
#' @return A phyloseq object with TSS normalized counts.
#' @references Bolstad et al., 2003.
norm.tss <- function(object) {
  otu <- otu_table(object)
  size <- colSums(otu)
  otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
  otu_table(object) <- otu_table(
    otu_normed,
    taxa_are_rows = taxa_are_rows(object)
  )
  object
}

#' CSS Normalization (Cumulative Sum Scaling)
#' @description Assumes equivalent count distributions for low abundant genes across samples up to a threshold. Scales the invariant segment of each sample’s count distribution.
#' @param ps A phyloseq object.
#' @param feature_category A character string specifying the feature category ("all", "iqlr", "zero", "lvha").
#' @return A phyloseq object with CSS normalized counts.
#' @references Paulson et al., 2013.
norm.css <- function(ps, feature_category = c("all", "iqlr", "zero", "lvha")) {
  ps <- remove_zero_negative_count_samples(ps)
  ps <- process_data_with_feature_category(ps, feature_category)
  
  if (nrow(otu_table(ps)) == 0 || ncol(otu_table(ps)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  
  cat("Initial OTU table dimensions: ", dim(otu_table(ps)), "\n")
  raw <- as(otu_table(ps), "matrix")
  raw <- raw + 1  # Add pseudocount
  dat.DGE <- DGEList(counts = raw)
  dat.DGE <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- dat.DGE$samples$norm.factors
  dat.normed <- t(t(raw) / scaling.factor)
  otu_table(ps) <- otu_table(dat.normed, taxa_are_rows = TRUE)
  set_nf(ps, scaling.factor)
}

#' RLE Normalization (Relative Log Expression)
#' @description Assumes most features are not differential and uses relative abundances to calculate the normalization factor.
#' @param object A phyloseq or otu_table object.
#' @param locfunc A function to compute the location statistic (default is median).
#' @param type A character string specifying the type of normalization ("poscounts" or "ratio").
#' @param geo_means A vector of geometric means for each feature.
#' @param control_genes A vector of control genes.
#' @return A phyloseq or otu_table object with normalization factors.
#' @references Anders and Huber, 2010.
norm.rle <- function(object, locfunc = stats::median, type = c("poscounts", "ratio"), geo_means = NULL, control_genes = NULL) {
  stopifnot(class(object) %in% c("phyloseq", "otu_table"))
  type <- match.arg(type, c("poscounts", "ratio"))
  
  geo_means <- ifelse(is.null(geo_means), substitute(), geo_means)
  control_genes <- ifelse(is.null(control_genes), substitute(), control_genes)
  
  otu <- as(otu_table(object), "matrix")
  nf <- estimateSizeFactorsForMatrix(
    otu,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )
  object_nf <- set_nf(object, nf)
  
  object_nf
}

#' Set normalization factors
#' @param object A phyloseq or otu_table object.
#' @param nf A vector of normalization factors.
#' @return The object with normalization factors set.
set_nf <- function(object, nf) {
  names(nf) <- NULL
  if (inherits(object, "phyloseq")) {
    sample_data(object) <- cbind(sample_data(object), norm_factor = nf)
    ot <- otu_table(object)
    attr(ot, "norm_factor") <- nf
    otu_table(object) <- ot
  } else if (inherits(object, "otu_table")) {
    attr(object, "norm_factor") <- nf
  } else {
    stop("object must be a `phyloseq` or `otu_table` object")
  }
  object
}

# Example usage for TC normalization
# Assuming ps is your phyloseq object and sample_data(ps)$Animal.type contains your group labels
#ps = subset_samples(physeq, !is.na(Animal.type))
# result_TC <- normalization_set(ps, method = "TC", groups = sample_data(ps)$Animal.type)
# normalized_ps_TC <- result_TC$dat.normed
# scaling_factors_TC <- result_TC$scaling.factor

# Example for UQ normalization
# result_UQ <- normalization_set(ps, method = "UQ", groups = sample_data(ps)$Animal.type)
# normalized_ps_UQ <- result_UQ$dat.normed
# scaling_factors_UQ <- result_UQ$scaling.factor

# Example for Median normalization
# result_med <- normalization_set(ps, method = "med", groups = sample_data(ps)$Animal.type)
# normalized_ps_med <- result_med$dat.normed
# scaling_factors_med <- result_med$scaling.factor

# Example for DESeq normalization
# result_DESeq <- normalization_set(ps, method = "DESeq", groups = sample_data(ps)$Animal.type)
# normalized_ps_DESeq <- result_DESeq$dat.normed
# scaling_factors_DESeq <- result_DESeq$scaling.factor

# Example for PoissonSeq normalization
# result_PoissonSeq <- normalization_set(ps, method = "PoissonSeq")
# normalized_ps_PoissonSeq <- result_PoissonSeq$dat.normed
# scaling_factors_PoissonSeq <- result_PoissonSeq$scaling.factor

# Example for Quantile normalization
# result_QN <- normalization_set(ps, method = "QN")
# normalized_ps_QN <- result_QN$dat.normed
# scaling_factors_QN <- result_QN$scaling.factor

# Example for SVA normalization
# result_SVA <- normalization_set(ps, method = "SVA", groups = sample_data(ps)$Animal.type)
# normalized_ps_SVA <- result_SVA$dat.normed
# scaling_factors_SVA <- result_SVA$scaling.factor

# Example for RUVg normalization
# result_RUVg <- normalization_set(ps, method = "RUVg", groups = sample_data(ps)$Animal.type)
# normalized_ps_RUVg <- result_RUVg$dat.normed
# scaling_factors_RUVg <- result_RUVg$scaling.factor

# Example for RUVs normalization
# result_RUVs <- normalization_set(ps, method = "RUVs", groups = sample_data(ps)$Animal.type)
# normalized_ps_RUVs <- result_RUVs$dat.normed
# scaling_factors_RUVs <- result_RUVs$scaling.factor
# ot<-normalized_ps_RUVs@otu_table
# write.csv(ot,"ot.csv")

# Example for RUVr normalization
# result_RUVr <- normalization_set(ps, method = "RUVr", groups = sample_data(ps)$Animal.type)
# normalized_ps_RUVr <- result_RUVr$dat.normed
# scaling_factors_RUVr <- result_RUVr$scaling.factor

# Example for TMM normalization
# result_TMM <- normalization_set(ps, method = "TMM", groups = sample_data(ps)$Animal.type)
# normalized_ps_TMM <- result_TMM$dat.normed
# scaling_factors_TMM <- result_TMM$scaling.factor

# Example for CLR normalization
# result_clr <- normalization_set(ps, method = "clr")
# normalized_ps_clr <- result_clr$dat.normed
# scaling_factors_clr <- result_clr$scaling.factor

# Example for Rarefying
# result_rar <- normalization_set(ps, method = "rar")
# normalized_ps_rar <- result_rar$dat.normed
# scaling_factors_rar <- result_rar$scaling.factor

# Example for CSS normalization
# result_css <- normalization_set(ps, method = "css")
# normalized_ps_css <- result_css$dat.normed
# scaling_factors_css <- result_css$scaling.factor

# Example for TSS normalization
# result_tss <- normalization_set(ps, method = "tss")
# normalized_ps_tss <- result_tss$dat.normed
# scaling_factors_tss <- result_tss$scaling.factor

# Example for RLE normalization
# result_rle <- normalization_set(ps, method = "rle")
# normalized_ps_rle <- result_rle$dat.normed
# scaling_factors_rle <- result_rle$scaling.factor

# Save the scaling factors to a file (example for DESeq normalization)
# write.csv(scaling_factors_DESeq, file = "scaling_factors_DESeq.csv", row.names = FALSE)
