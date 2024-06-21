#' Run Variance Stabilizing Transformation (VST) Analysis
#'
#' This function performs a variance stabilizing transformation (VST) on the counts data of a phyloseq object
#' using DESeq2, and returns a phyloseq object with the transformed data.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param design_formula A formula specifying the design of the experiment. Default is ~ sample_name.
#' @param output_file A character string specifying the output file name for the transformed DESeqDataSet object. Default is NULL.
#' @param pseudocount A numeric value specifying the pseudocount to add to the counts data before transformation. Default is 1.
#' @return A phyloseq object with the VST-transformed counts data.
#' @examples
#' # Perform VST analysis on a phyloseq object
#' transformed_physeq <- run_vst_analysis(physeq_ITS)
#' @export
run_vst_analysis <- function(physeq, design_formula = ~ sample_name, output_file = NULL, pseudocount = 1) {
  # Load necessary libraries
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
  deseq_data <- phyloseq_to_deseq2(physeq, design_formula)
  
  # Ensure all variables in the design formula are factors
  design_vars <- all.vars(design_formula)
  for (var in design_vars) {
    colData(deseq_data)[[var]] <- factor(colData(deseq_data)[[var]])
  }
  
  # Ensure counts are integers and add a small pseudocount
  counts_matrix <- counts(deseq_data)
  counts_matrix <- round(counts_matrix + pseudocount)
  
  # Check if there are any non-finite values in the counts matrix
  if (any(!is.finite(counts_matrix))) {
    stop("Counts matrix contains non-finite values.")
  }
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                colData = colData(deseq_data),
                                design = design_formula)
  
  # Perform variance stabilizing transformation (VST)
  vst_dds <- varianceStabilizingTransformation(dds)
  
  # Save the result if an output file is specified
  if (!is.null(output_file)) {
    saveRDS(vst_dds, file = output_file)
    cat("Transformed DESeqDataSet object saved to:", output_file, "\n")
  }
  
  # Create a phyloseq object with the transformed OTU table
  vst_physeq <- phyloseq(phyloseq::otu_table(assay(vst_dds), taxa_are_rows = TRUE),
                         phyloseq::sample_data(physeq),
                         phyloseq::tax_table(physeq))
  
  return(vst_physeq)
}

# Example usage:
# Perform VST analysis on a phyloseq object
# transformed_physeq <- run_vst_analysis(physeq_ITS)
