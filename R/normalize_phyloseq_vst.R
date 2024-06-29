#' Variance Stabilizing Transformation (VST)
#'
#' This function performs VST normalization on a phyloseq object using the DESeq2 package.
#' It handles zero, negative, and NA values by removing them and informing the user of the modifications.
#'
#' @param physeq A phyloseq object.
#' @param design_formula A formula specifying the design for DESeq2.
#' @param feature_category Character string specifying which features to use as the divisor for the geometric mean calculation.
#' @param min_counts Minimum number of counts required for a sample to be retained.
#' @param output_file An optional file path to save the transformed DESeqDataSet object.
#' @param pseudocount A numeric value to add to counts before transformation to avoid log(0).
#' @return A phyloseq object with VST normalization.
#' @export
#' @examples
#' data(GlobalPatterns)
#' ps_vst <- normalize_phyloseq_vst(GlobalPatterns, design_formula = ~ SampleType, feature_category = "iqlr", min_counts = 1000)
normalize_phyloseq_vst <- function(physeq, design_formula = ~ sample_name, feature_category = c("all", "iqlr", "zero", "lvha"), min_counts = 1, output_file = NULL, pseudocount = 1) {
  physeq <- remove_zero_negative_count_samples(physeq) # Remove samples with zero or negative counts
  physeq <- process_data_with_feature_category(physeq, feature_category) # Process data with the specified feature category
  
  # Check if the OTU table has valid dimensions after processing
  if (nrow(otu_table(physeq)) == 0 || ncol(otu_table(physeq)) == 0) {
    stop("OTU table has non-zero dimensions after processing. Ensure the data contains valid counts.")
  }
  
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required but not installed.")
  }
  library(phyloseq)
  library(DESeq2)
  
  # Check if physeq is a valid phyloseq object
  if (!is(physeq, "phyloseq")) {
    stop("Input 'physeq' must be a valid phyloseq object.")
  }
  
  # Add sample_name to the sample data if not present
  if (!"sample_name" %in% colnames(sample_data(physeq))) {
    sample_data(physeq)$sample_name <- sample_names(physeq)
  }
  
  # Convert phyloseq to DESeqDataSet
  deseq_data <- phyloseq::phyloseq_to_deseq2(physeq, design_formula)
  
  # Ensure all variables in the design formula are factors
  design_vars <- all.vars(design_formula)
  for (var in design_vars) {
    colData(deseq_data)[[var]] <- factor(colData(deseq_data)[[var]])
  }
  
  # Ensure counts are integers and add a small pseudocount
  counts_matrix <- DESeq2::counts(deseq_data)
  counts_matrix <- round(counts_matrix + pseudocount)
  
  # Check for non-finite values and replace them with 0
  counts_matrix[!is.finite(counts_matrix)] <- 0
  
  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_matrix,
                                        colData = colData(deseq_data),
                                        design = design_formula)
  
  # Estimate dispersions using gene-wise estimates as final estimates
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  
  # Perform variance stabilizing transformation (VST)
  vst_dds <- DESeq2::varianceStabilizingTransformation(dds)
  
  # Save the result if an output file is specified
  if (!is.null(output_file)) {
    saveRDS(vst_dds, file = output_file)
    cat("Transformed DESeqDataSet object saved to:", output_file, "\n")
  }
  
  # Create a phyloseq object with the transformed OTU table
  vst_physeq <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(assay(vst_dds)), taxa_are_rows = TRUE),
                                   phyloseq::sample_data(physeq),
                                   phyloseq::tax_table(physeq))
  
  return(vst_physeq)
}

# Example;
# physeq <- phy
# physeq_vst <- normalize_phyloseq_vst(physeq, design_formula = ~treatment, feature_category = "iqlr", min_counts = 1000)
