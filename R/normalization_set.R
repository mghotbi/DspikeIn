#' Overview of Normalization Methods
#'
#' @description This documentation provides an overview of normalization methods used for microbiome analysis.
#'
#' - **For Differential Abundance Analysis**: `DESeq`, `TMM`, or `CSS` normalization methods handle compositional biases and library size differences.
#' - **For Compositional Data**: `CLR` normalization accounts for compositional structure by transforming the data into log-ratio format.
#' - **For Simplicity and Ease of Use**: `TC`, `UQ`, or `Median` normalization methods are quick but may not be as robust.
#'
#' @section Normalization Use Cases:
#' These normalization methods are commonly used in microbiome analysis to ensure fair comparisons across samples.
#'
#' @keywords normalization microbiome NULL
#' @title Calculate Geometric Mean
#' 
#' @name gm_mean
#'
#' @description This function calculates the geometric mean of a numeric vector.
#' It removes non-positive and NA values by default.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values (NAs) be removed? Defaults to TRUE.
#' @return Geometric mean of x, or NA if no valid values are present.
#' @examples
#' vec <- c(1, 10, 100, 1000)
#' gm_mean(vec)
#' @export
gm_mean <- function(x, na.rm = TRUE) {
  valid_x <- x[x > 0 & !is.na(x)]
  if (length(valid_x) == 0) {
    return(NA)
  }
  exp(sum(log(valid_x), na.rm = na.rm) / length(valid_x))
}
# ------------------------------------------------------------------------
#' @title Set Normalization Factors in the Sample Data of the Phyloseq Object
#' @name set_nf
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @param scaling.factor A vector of normalization factors.
#' @return A phyloseq object with updated sample data.
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   data("physeq_ITSOTU", package = "DspikeIn")
#'
#'   # Create normalization factors (e.g., all ones)
#'   nf <- rep(1, phyloseq::nsamples(physeq_ITSOTU))
#'
#'   # Apply normalization factors
#'   physeq_ITSOTU <- DspikeIn::set_nf(physeq_ITSOTU, scaling.factor = nf)
#'
#'   # Check the updated sample data
#'   head(phyloseq::sample_data(physeq_ITSOTU)$norm_factors)
#' }
#' @export
set_nf <- function(obj, scaling.factor) {
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  }

  if (length(scaling.factor) != ncol(get_otu_table(obj))) {
    stop("Length of scaling.factor must match the number of samples.")
  }

  if (inherits(obj, "phyloseq")) {
    phyloseq::sample_data(obj)$norm_factors <- scaling.factor
  } else {
    SummarizedExperiment::colData(obj)$norm_factors <- scaling.factor
  }

  return(obj)
}

# -----------------------------------------------------------
#' Tidy a Phyloseq Object and Remove Zero/Negative Count Samples
#'
#' @title Tidy a Phyloseq or TreeSummarizedExperiment Object
#' @name tidy_phyloseq_tse
#' @description Cleans and standardizes a microbiome dataset, supporting both
#' `phyloseq` and `TreeSummarizedExperiment`. Performs:
#' - Standardization of taxonomic ranks (if available)
#' - Removal of leading/trailing whitespace in taxa names
#' - Filtering out zero-count taxa
#' - Exclusion of "Chloroplast" and "Mitochondria" classifications (if applicable)
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A cleaned and tidied object of the same class.
#'
#' @details This function standardizes taxonomic ranks, removes unnecessary whitespace,
#' and filters unwanted classifications, ensuring consistency for downstream analysis.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   tidy_physeq <- tidy_phyloseq_tse(physeq_16SOTU)
#' }
#'
#' @importFrom phyloseq prune_taxa
#' @importFrom SummarizedExperiment assay
tidy_phyloseq_tse <- function(obj) {
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  # Extract taxonomy table
  tax_data <- get_tax_table(obj)

  if (is.null(tax_data) || ncol(tax_data) == 0) {
    warning("No taxonomy table found. Returning unmodified object.")
    return(obj)
  }

  # Fix taxa names by removing "__" prefixes (e.g., "k__Bacteria" becomes "Bacteria")
  tax_data <- as.data.frame(lapply(tax_data, function(col) gsub("[a-z]__\\s*", "", col)), stringsAsFactors = FALSE)

  # Define standard taxonomic ranks
  required_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  available_ranks <- intersect(required_ranks, colnames(tax_data))

  # Ensure taxonomic table maintains structure
  if (length(available_ranks) > 0) {
    tax_data <- tax_data[, available_ranks, drop = FALSE]
  } else {
    warning("No standard taxonomic ranks found. Proceeding with available data.")
  }

  # Trim whitespace from taxonomic names
  tax_data <- as.data.frame(lapply(tax_data, trimws), stringsAsFactors = FALSE)

  # Extract OTU table
  otu_matrix <- get_otu_table(obj)

  if (is.null(otu_matrix) || nrow(otu_matrix) == 0) {
    warning("No OTU table found. Returning unmodified object.")
    return(obj)
  }

  # Identify and remove zero-count taxa
  non_zero_taxa <- rowSums(otu_matrix) > 0

  # Remove unwanted taxa based on taxonomy
  if ("Class" %in% colnames(tax_data)) {
    non_zero_taxa <- non_zero_taxa & tax_data$Class != "Chloroplast"
  }
  if ("Family" %in% colnames(tax_data)) {
    non_zero_taxa <- non_zero_taxa & tax_data$Family != "Mitochondria"
  }

  # Apply filtering to the object
  if (inherits(obj, "phyloseq")) {
    obj <- phyloseq::prune_taxa(rownames(otu_matrix)[non_zero_taxa], obj)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- obj[rownames(otu_matrix)[non_zero_taxa], ]
  }

  return(obj)
}

#' @title Remove Samples with Zero, Negative Counts, or NA Values and Add Pseudocount
#' @name  remove_zero_negative_count_samples
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @param pseudocount A numeric value to add to avoid zero counts.
#' @return A phyloseq object with filtered and adjusted OTU table.
#' @importFrom phyloseq otu_table prune_samples prune_taxa sample_sums
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   library(DspikeIn)
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Remove samples with zero/negative/NA counts and add pseudocount
#'   cleaned_ps <- remove_zero_negative_count_samples(
#'     physeq_16SOTU,
#'     pseudocount = 1e-6
#'   )
#' }
#' @export
remove_zero_negative_count_samples <- function(obj, pseudocount = 1e-6) {
  if (inherits(obj, "phyloseq")) {
    otu <- as(phyloseq::otu_table(obj), "matrix")

    # Identify samples with zero, negative counts, or NA values
    zero_negative_count_samples <- phyloseq::sample_sums(obj) <= 0
    na_count_samples <- apply(otu, 2, function(x) any(is.na(x)))
    samples_to_remove <- zero_negative_count_samples | na_count_samples

    # Prune the samples that meet the removal criteria
    if (any(samples_to_remove)) {
      message(sprintf("Removing %d samples with zero, negative counts, or NA values.", sum(samples_to_remove)))
      obj <- phyloseq::prune_samples(!samples_to_remove, obj)
      otu <- as(phyloseq::otu_table(obj), "matrix")
    }

    # Remove features with zero counts across all samples
    zero_rows <- rowSums(otu) == 0
    if (any(zero_rows)) {
      message(sprintf("Removing %d features with zero counts across all samples.", sum(zero_rows)))
      obj <- phyloseq::prune_taxa(!zero_rows, obj)
    }

    # Add pseudocount and return
    otu <- as(phyloseq::otu_table(obj), "matrix") + pseudocount
    phyloseq::otu_table(obj) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
    return(obj)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    #  Extract count data
    otu <- SummarizedExperiment::assay(obj)

    # Identify and remove zero-count taxa
    zero_rows <- rowSums(otu) == 0
    if (any(zero_rows)) {
      message(sprintf("Removing %d features with zero counts across all samples.", sum(zero_rows)))
      # Ensure both `assay()` and `rowData()` remain synchronized
      obj <- obj[!zero_rows, ] # Removes both count data and metadata
    }

    # Add pseudocount
    otu <- SummarizedExperiment::assay(obj) + pseudocount
    SummarizedExperiment::assay(obj) <- otu

    return(obj)
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}


# -----------------------------------------------------------
#' @title Convert Categorical Columns to Factors in Sample Data
#' @name convert_categorical_to_factors
#' @importFrom phyloseq sample_data
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing microbial data.
#' @return A phyloseq object with updated sample data.
#' @examples
#' data("physeq_16SOTU", package = "DspikeIn")
#' ps_factor <- convert_categorical_to_factors(physeq_16SOTU)
#' @export
convert_categorical_to_factors <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    sample_data_df <- as(phyloseq::sample_data(obj), "data.frame")
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    sample_data_df <- as.data.frame(SummarizedExperiment::colData(obj))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }

  for (col in colnames(sample_data_df)) {
    if (is.character(sample_data_df[[col]]) || is.factor(sample_data_df[[col]])) {
      sample_data_df[[col]] <- as.factor(sample_data_df[[col]])
    }
  }

  if (inherits(obj, "phyloseq")) {
    phyloseq::sample_data(obj) <- phyloseq::sample_data(sample_data_df)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    SummarizedExperiment::colData(obj) <- sample_data_df
  }

  return(obj)
}

# -----------------------------------------------------------
#' @title Create a List from a Phyloseq or TSE Object
#' @name create_list
#' @importFrom phyloseq otu_table sample_data tax_table
#' @param obj A phyloseq or TreeSummarizedExperiment object.
#' @return A list containing the DGE list and updated phyloseq object.
create_list <- function(obj) {
  # Convert if needed
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
    if (!inherits(obj, "phyloseq")) stop("Failed to convert TSE to phyloseq.")
  }

  # Ensure it's a valid phyloseq object
  if (!inherits(obj, "phyloseq")) stop("Input must be a phyloseq object.")

  # Extract components safely
  counts <- get_otu_table(obj)
  sample_data_df <- get_sample_data(obj)

  if (is.null(counts) || nrow(counts) == 0) stop("OTU table extraction failed.")
  if (is.null(sample_data_df)) stop("Sample data extraction failed.")

  # Create lib size and normalization factors
  lib.size <- colSums(counts)
  norm.factors <- rep(1, ncol(counts))

  sample_data_df$lib.size <- lib.size
  sample_data_df$norm.factors <- norm.factors

  # Ensure correct structure
  dge_base <- list(
    counts = counts,
    samples = sample_data_df,
    group = sample_data_df$group,
    lib.size = lib.size,
    norm.factors = norm.factors
  )

  return(list(dge_list = dge_base, phyloseq_obj = obj))
}


# -----------------------------------------------------------
#' @title Apply the Selected Normalization Method to the Phyloseq and TSE Objects
#' @name normalization_set
#' @param obj A phyloseq object.
#' @param method A character string specifying the normalization method ("TC", "UQ", "med", "DESeq", "Poisson", "QN", "TMM", "clr", "rar", "css", "tss", "rle").
#' @param groups A column name of group labels from sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
#' @examples
#' # Example with a phyloseq object
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   ps <- physeq_16SOTU
#'   result_phyloseq <- normalization_set(ps, method = "TC", groups = "Host.species")
#'   head(result_phyloseq$scaling.factor)
#'   normed_physeq <- result_phyloseq$dat.normed
#' }
#'
#' # Example with a TreeSummarizedExperiment (TSE) object
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   result_tse <- normalization_set(tse_16SOTU, method = "clr")
#'   head(result_tse$scaling.factor)
#'   normed_tse <- result_tse$dat.normed
#' }
#'
#' # For a full comparison of all normalization methods, see the vignette:
#' # vignette("DspikeIn-with-Phyloseq", package = "DspikeIn")
#' # vignette("DspikeIn-with-TSE", package = "DspikeIn")
#'
#' @export
normalization_set <- function(obj, method, groups = NULL) {
  # Convert TSE to phyloseq if needed
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  if (phyloseq::nsamples(obj) == 0) stop("The phyloseq object contains no samples.")
  if (length(phyloseq::sample_names(obj)) == 0) stop("Sample names are missing in the phyloseq object.")

  obj <- remove_zero_negative_count_samples(obj)
  obj <- convert_categorical_to_factors(obj)

  # Ensure groups is a single column name from sample data
  if (!is.null(groups)) {
    if (!is.character(groups) || length(groups) != 1 || !groups %in% colnames(phyloseq::sample_data(obj))) {
      stop("'groups' should be a single column name from sample data.")
    }
  }

  result <- switch(method,
    "TC" = norm.TC(obj, groups),
    "UQ" = norm.UQ(obj, groups),
    "med" = norm.med(obj, groups),
    "DESeq" = norm.DESeq(obj, groups),
    "Poisson" = norm.Poisson(obj, groups),
    "QN" = norm.QN(obj),
    "TMM" = norm.TMM(obj, groups),
    "clr" = norm.clr(obj),
    "rar" = norm.rar(obj),
    "css" = norm.css(obj),
    "tss" = norm.tss(obj),
    "rle" = norm.rle(obj),
    stop("Invalid normalization method")
  )

  dat.normed <- result$dat.normed
  scaling.factor <- result$scaling.factor

  return(list(dat.normed = dat.normed, scaling.factor = scaling.factor))
}
# -----------------------------------------------------------
#' @title TC Normalization (Total Count Scaling)
#' @name norm.TC
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom edgeR DGEList
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param groups A string specifying the grouping variable in sample data.
#'
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.TC <- function(obj, groups) {
  # Convert TSE to phyloseq if necessary
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  # Ensure it's a phyloseq object
  if (!inherits(obj, "phyloseq")) {
    stop("Input must be a phyloseq or TreeSummarizedExperiment object.")
  }

  # Create list from phyloseq
  physeq_list <- create_list(obj)
  obj <- physeq_list$phyloseq_obj
  dat.DGE <- physeq_list$dge_list

  # Compute scaling factor
  scaling.factor <- dat.DGE$samples$lib.size / 1e6

  # Apply total count normalization
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)

  # Update OTU table
  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)

  # Set normalization factors in sample data
  obj <- set_nf(obj, scaling.factor)

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}


# -----------------------------------------------------------
#' @title UQ Normalization (Upper Quartile)
#' @name norm.UQ
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom stats quantile
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.UQ <- function(obj, groups) {
  # Create a custom list
  physeq_list <- create_list(obj)
  obj <- physeq_list$phyloseq_obj
  dat.DGE <- physeq_list$dge_list

  # Calculate the upper quartile factor
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))

  # Handle cases where q.factor might be zero
  if (any(q.factor == 0)) stop("One or more upper quartile factors are zero.")

  scaling.factor <- q.factor / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)

  # Update the OTU table in the phyloseq object
  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)

  # Set normalization factors
  obj <- set_nf(obj, scaling.factor)

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title Median Normalization
#' @name norm.med
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom stats median
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.med <- function(obj, groups) {
  obj <- create_list(obj)$phyloseq_obj
  dat.DGE <- create_list(obj)$dge_list
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6
  dat.normed <- t(t(dat.DGE$counts) / scaling.factor)
  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  obj <- set_nf(obj, scaling.factor)
  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' DESeq Normalization
#'
# -----------------------------------------------------------
#' @title DESeq Normalization with Pseudocount and Integer Conversion
#' @name norm.DESeq
#' @importFrom phyloseq otu_table taxa_are_rows sample_data
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq sizeFactors counts estimateSizeFactors
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @param groups A string specifying the grouping variable in sample data.
#' @param pseudocount A numeric value added to avoid zeros in the dataset.
#' @return A list containing the normalized object (same format as input) and scaling factors.
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Example 1: phyloseq input (subset to Animal.type == "Frog")
#'   physeq_frog <- phyloseq::subset_samples(physeq_16SOTU, Animal.type == "Frog")
#'   result_DESeq_phy <- norm.DESeq(physeq_frog, groups = "Animal.type", pseudocount = 1)
#'
#'   # Example 2: TSE input (convert and subset to Animal.type == "Frog")
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   col_meta <- SummarizedExperiment::colData(tse_16SOTU)
#'   tse_frog <- tse_16SOTU[, which(col_meta$Animal.type == "Frog")]
#'   result_DESeq_tse <- norm.DESeq(tse_frog, groups = "Animal.type", pseudocount = 1)
#' }
#' @export
norm.DESeq <- function(obj, groups, pseudocount = 1) {
  # Ensure input is valid
  if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
    stop("Error: Input must be a 'phyloseq' or 'TreeSummarizedExperiment' object.")
  }

  # Patch for broken DESeq2 S4 class hierarchy in some R installs
  if (!methods::isClass("ExpData")) {
    setClass("ExpData", contains = "VIRTUAL")
  }

  # Convert TSE to phyloseq if needed
  input_class <- class(obj)
  if (inherits(obj, "TreeSummarizedExperiment")) {
    obj <- convert_tse_to_phyloseq(obj)
  }

  # Clean Data
  obj <- remove_zero_negative_count_samples(obj)
  obj <- convert_categorical_to_factors(obj)

  # Extract Count Table
  raw <- as(phyloseq::otu_table(obj), "matrix") + pseudocount
  raw <- round(raw) # required by DESeq2

  # Extract Sample Data
  sample_data_df <- as.data.frame(phyloseq::sample_data(obj))
  if (!groups %in% colnames(sample_data_df)) {
    stop("Error: Specified 'groups' column not found in sample data.")
  }
  sample_data_df[[groups]] <- as.factor(sample_data_df[[groups]])

  # Create DESeq2 design formula dynamically
  condition <- data.frame(Condition = sample_data_df[[groups]])
  rownames(condition) <- colnames(raw)

  if (nlevels(condition$Condition) == 1) {
    warning("Only one group detected! Using design ~ 1.")
    design_formula <- ~1
  } else {
    design_formula <- ~Condition
  }

  # Create DESeq2 object
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw,
    colData = condition,
    design = design_formula
  )

  # Run DESeq normalization
  if (nlevels(condition$Condition) == 1) {
    dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  } else {
    dat.DGE <- DESeq2::DESeq(dat.DGE, fitType = "local")
  }

  # Extract results
  scaling.factor <- DESeq2::sizeFactors(dat.DGE)
  dat.normed <- DESeq2::counts(dat.DGE, normalized = TRUE)

  # Update phyloseq object
  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  obj <- set_nf(obj, scaling.factor)

  # Convert back to TSE if needed
  if (input_class == "TreeSummarizedExperiment") {
    obj <- convert_phyloseq_to_tse(obj)
  }

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title Quantile Normalization (QN) for phyloseq object
#' @name norm.QN
#' @importFrom phyloseq otu_table taxa_are_rows taxa_names tax_table
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param filter Logical, whether to filter low counts.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.QN <- function(obj, filter = FALSE) {
  otu <- as(phyloseq::otu_table(obj), "matrix")

  if (!is.numeric(otu)) {
    stop("OTU table must contain numeric values.")
  }

  taxa_names_original <- phyloseq::taxa_names(obj)
  tax_table_original <- phyloseq::tax_table(obj)

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
  phyloseq::otu_table(obj) <- phyloseq::otu_table(normalized, taxa_are_rows = TRUE)

  # for QN we don't need scaling factors
  scaling.factor <- NULL

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title Poisson Normalization and Differential Abundance Function
#' @name norm.Poisson
#' @importFrom phyloseq otu_table sample_data taxa_are_rows
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @importFrom stats model.matrix
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param group_var A string specifying the grouping variable in sample data (if phyloseq object).
#' @param pseudocount A numeric value added to avoid division by zero.
#' @return A list containing the normalized data, scaling factor, and differential abundance results.
norm.Poisson <- function(obj, group_var = NULL, pseudocount = 1e-6) {
  obj <- remove_zero_negative_count_samples(obj)

  prepare_data <- function(obj, group_var) {
    if (inherits(obj, "phyloseq")) {
      raw_otu <- as(phyloseq::otu_table(obj), "matrix")
      sample_data_df <- as(phyloseq::sample_data(obj), "data.frame")

      if (!is.null(group_var)) {
        y <- as.numeric(as.factor(sample_data_df[[group_var]]))
      } else {
        y <- rep(1, ncol(raw_otu))
      }
    } else if (is.matrix(obj)) {
      raw_otu <- obj
      y <- rep(1, ncol(raw_otu))
    } else {
      stop("Input obj must be a phyloseq object or a matrix.")
    }
    list(raw_otu = raw_otu, sample_data_df = sample_data_df, y = y)
  }

  prepared_data <- prepare_data(obj, group_var)
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
  design <- stats::model.matrix(~y)
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit)
  topTags <- edgeR::topTags(lrt, n = nrow(raw_otu))

  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  obj <- set_nf(obj, scaling.factor)
  return(list(dat.normed = obj, scaling.factor = scaling.factor, differential_abundance = topTags))
}

# -----------------------------------------------------------
#' @title TMM Normalization (Trimmed Mean of M component)
#' @name norm.TMM
#' @importFrom phyloseq otu_table taxa_are_rows sample_data
#' @importFrom edgeR DGEList calcNormFactors
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param groups A string specifying the grouping variable in sample data.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.TMM <- function(obj, groups) {
  otu_table_matrix <- as(phyloseq::otu_table(obj), "matrix")

  zero_rows <- rowSums(otu_table_matrix) == 0
  if (any(zero_rows)) {
    otu_table_matrix <- otu_table_matrix[!zero_rows, ]
  }

  zero_cols <- colSums(otu_table_matrix) == 0
  if (any(zero_cols)) {
    otu_table_matrix <- otu_table_matrix[, !zero_cols]
  }

  sample_data_df <- as(phyloseq::sample_data(obj), "data.frame")
  group <- sample_data_df[[groups]]
  dge <- edgeR::DGEList(counts = otu_table_matrix, group = group)

  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  scaling.factor <- dge$samples$norm.factors
  dat.normed <- t(t(dge$counts) / scaling.factor)

  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = phyloseq::taxa_are_rows(obj))
  phyloseq::sample_data(obj)$norm_factors <- scaling.factor

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title CLR Normalization (Centered Log-Ratio Transformation)
#' @name norm.clr
#' @importFrom phyloseq transform_sample_counts nsamples
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.clr <- function(obj) {
  obj <- remove_zero_negative_count_samples(obj)

  gm_mean <- function(x) exp(sum(log(x[x > 0])) / length(x))

  obj_clr <- phyloseq::transform_sample_counts(obj, function(x) log(x / gm_mean(x)))
  scaling.factor <- rep(1, phyloseq::nsamples(obj_clr))

  return(list(dat.normed = obj_clr, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title Rarefying
#' @name norm.rar
#' @importFrom phyloseq rarefy_even_depth sample_sums
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.rar <- function(obj) {
  obj <- remove_zero_negative_count_samples(obj)

  obj_rarefied <- phyloseq::rarefy_even_depth(obj, rngseed = 123)
  scaling.factor <- phyloseq::sample_sums(obj_rarefied)

  return(list(dat.normed = obj_rarefied, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title TSS Normalization (Total Sum Scaling)
#' @name norm.tss
#' @importFrom phyloseq otu_table taxa_are_rows nsamples
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @return A list containing the normalized phyloseq object and scaling factor
norm.tss <- function(obj) {
  obj <- remove_zero_negative_count_samples(obj)

  otu <- phyloseq::otu_table(obj)
  size <- colSums(otu)
  otu_normed <- sweep(otu, MARGIN = 2, STATS = size, FUN = "/")
  phyloseq::otu_table(obj) <- phyloseq::otu_table(otu_normed, taxa_are_rows = phyloseq::taxa_are_rows(obj))

  scaling.factor <- rep(1, phyloseq::nsamples(obj))

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title CSS Normalization (Cumulative Sum Scaling)
#' @name norm.css
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom edgeR DGEList calcNormFactors
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.css <- function(obj) {
  obj <- remove_zero_negative_count_samples(obj)

  raw <- as(phyloseq::otu_table(obj), "matrix") + 1e-6
  dge <- edgeR::DGEList(counts = raw)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  scaling.factor <- dge$samples$norm.factors
  dat.normed <- t(t(raw) / scaling.factor)

  phyloseq::otu_table(obj) <- phyloseq::otu_table(dat.normed, taxa_are_rows = TRUE)
  obj <- set_nf(obj, scaling.factor)

  return(list(dat.normed = obj, scaling.factor = scaling.factor))
}

# -----------------------------------------------------------
#' @title RLE Normalization (Relative Log Expression)
#' @name norm.rle
#' @importFrom phyloseq otu_table taxa_are_rows
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom stats median
#'
#' @param obj A Phyloseq or TreeSummarizedExperiment objects.
#' @param locfunc A function to compute the location statistic (default is median).
#' @param type A character string specifying the type of normalization ("poscounts" or "ratio").
#' @param geo_means A vector of geometric means for each feature.
#' @param control_genes A vector of control genes.
#' @return A list containing the normalized phyloseq object and scaling factors.
norm.rle <- function(obj, locfunc = stats::median, type = c("poscounts", "ratio"), geo_means = NULL, control_genes = NULL) {
  type <- match.arg(type, c("poscounts", "ratio"))

  #  Extract OTU table correctly
  otu <- as(phyloseq::otu_table(obj), "matrix")

  #  Ensure geo_means is valid
  if (is.null(geo_means)) {
    geo_means <- apply(otu, 1, gm_mean)
  }

  #  Ensure control_genes is correctly formatted
  if (!is.null(control_genes)) {
    if (!is.numeric(control_genes) && !is.logical(control_genes)) {
      stop("Error: controlGenes must be a numeric or logical vector.")
    }

    if (is.numeric(control_genes)) {
      if (any(control_genes < 1 | control_genes > nrow(otu))) {
        stop("Error: controlGenes contains indices out of range.")
      }
    }
  } else {
    # Default: Use all genes as control genes
    control_genes <- rep(TRUE, nrow(otu)) # A logical vector selecting all genes
  }

  #  Call DESeq2 normalization
  nf <- DESeq2::estimateSizeFactorsForMatrix(
    otu,
    locfunc = locfunc,
    geoMeans = geo_means,
    controlGenes = control_genes,
    type = type
  )

  #  Update phyloseq object
  phyloseq::otu_table(obj) <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
  obj <- set_nf(obj, nf)

  return(list(dat.normed = obj, scaling.factor = nf))
}

# -----------------------------------------------------------
# Example usage for TC normalization
# ps is phyloseq object and sample_data(ps)$Animal.type
# contains your group labels
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
# result_med <- normalization_set(ps,
# method = "med", groups = "Host.species")
# normalized_ps_med <- result_med$dat.normed
# scaling_factors_med <- result_med$scaling.factor

# -----------------------------------------------------------
# Example for DESeq normalization
# ps_n <-remove_zero_negative_count_samples(ps)
# result_DESeq <- normalization_set(ps,
# method = "DESeq", groups = "Animal.type")
# normalized_ps_DESeq <- result_DESeq$dat.normed
# scaling_factors_DESeq <- result_DESeq$scaling.factor

# -----------------------------------------------------------
# Example for Poisson normalization
# result_Poisson <- normalization_set(ps,
# method = "Poisson", groups = "Host.genus")
# normalized_ps_Poisson <- result_Poisson$dat.normed
# scaling_factors_Poisson <- result_Poisson$scaling.factor

# -----------------------------------------------------------
# Example for Quantile normalization
# result_QN <- normalization_set(ps, method = "QN")
# normalized_ps_QN <- result_QN$dat.normed
# scaling_factors_QN <- result_QN$scaling.factor

# -----------------------------------------------------------
# Example for TMM normalization
# result_TMM <- normalization_set(ps,
# method = "TMM", groups = "Animal.type")
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
