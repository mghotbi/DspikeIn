#' Select ASVs Using Random Forest
#'
#' This function selects the most important ASVs (Amplicon Sequence Variants) based on a Random Forest model.
#' The selected ASVs are saved to a CSV file.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param response_var A character string specifying the response variable from the sample data.
#' @param minlib A numeric value specifying the minimum library size. Default is 15000.
#' @param prunescale A numeric value specifying the pruning scale for rare OTUs. Default is 0.0001.
#' @param ntree An integer specifying the number of trees to grow in the Random Forest. Default is 100.
#' @param n_top_predictors An integer specifying the number of top predictors to select. Default is 20.
#' @param output_csv A character string specifying the output CSV file name. Default is "randomforest_results.csv".
#' @param na_vars A character vector specifying the sample data variables to check for NA values. Default is NULL.
#' @return A pruned phyloseq object containing the selected ASVs.
#' @examples
#' # Ensure you have a phyloseq object `physeq`
#' # Use the name of the column in sample_data for response_var
#' rf_physeq <- RandomForest_selected_ASVs(physeq, response_var = "host.species", na_vars = c("Habitat", "Ecoregion_III", "Host_genus", "Animal_type"))
#' @export
RandomForest_selected_ASVs <- function(physeq, response_var, minlib = 15000, prunescale = 0.0001, ntree = 100, n_top_predictors = 20, output_csv = "randomforest_results.csv", na_vars = NULL) {
  
  # Load necessary libraries
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required but not installed.")
  }
  
  library(phyloseq)
  library(dplyr)
  library(randomForest)
  
  # Check if response_var exists in sample_data
  if (!response_var %in% colnames(sample_data(physeq))) {
    stop(paste("Response variable", response_var, "not found in sample data."))
  }
  
  # Ensure na_vars includes response_var
  if (!is.null(na_vars)) {
    na_vars <- unique(c(response_var, na_vars))
  } else {
    na_vars <- response_var
  }
  
  # Remove samples with NA in the specified variables
  sample_data_df <- data.frame(sample_data(physeq))
  for (var in na_vars) {
    sample_data_df <- sample_data_df[!is.na(sample_data_df[[var]]), ]
  }
  sample_data(physeq) <- sample_data(sample_data_df)
  
  # Convert categorical data to factors
  for (var in na_vars) {
    sample_data(physeq)[[var]] <- as.factor(sample_data(physeq)[[var]])
  }
  
  # Check if response variable has at least two classes
  if (length(unique(sample_data_df[[response_var]])) < 2) {
    stop("The response variable needs at least two classes for classification.")
  }
  
  # Remove taxa with any NA values in tax_table
  physeq <- subset_taxa(physeq, apply(tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
  
  # Remove samples with zero reads after subsetting
  physeq <- prune_samples(sample_sums(physeq) > 0, physeq)
  
  # Prune out rare OTUs by mean relative abundance set by prunescale
  tax.mean <- taxa_sums(physeq) / nsamples(physeq)
  s.prune <- prune_taxa(tax.mean > prunescale * minlib, physeq)
  
  # Replace empty or NA values in tax_table with "unidentified"
  tax_table(s.prune) <- replace(tax_table(s.prune), is.na(tax_table(s.prune)), "unidentified")
  tax_table(s.prune) <- replace(tax_table(s.prune), tax_table(s.prune) == "", "unidentified")
  
  predictors <- t(otu_table(s.prune))
  
  # Make one column for our outcome/response variable 
  response <- as.factor(sample_data(s.prune)[[response_var]])
  
  # Combine into one data frame
  rf.data <- data.frame(response, predictors)
  
  # Remove rows with NA values from rf.data
  rf.data <- rf.data[complete.cases(rf.data), ]
  
  # Check again if response variable has at least two classes after removing NAs
  if (length(unique(rf.data$response)) < 2) {
    stop("The response variable needs at least two classes for classification after removing NAs.")
  }
  
  # Random Forest model/for reproducibility
  set.seed(2)
  sal.classify <- randomForest(response ~ ., data = rf.data, ntree = ntree, keep.forest = FALSE, proximity = FALSE, importance = TRUE)
  
  # Extract variable importance
  imp <- importance(sal.classify)
  imp_df <- data.frame(predictors = rownames(imp), imp)
  imp_df <- arrange(imp_df, desc(MeanDecreaseGini))
  imp_df$predictors <- factor(imp_df$predictors, levels = imp_df$predictors)
  
  # Select the top n_top_predictors 
  imp_top <- imp_df[1:n_top_predictors, ]
  
  # Subset selected ASVs from phyloseq object
  otunames <- imp_top$predictors
  r <- rownames(tax_table(s.prune)) %in% otunames
  selected_phyloseq <- prune_taxa(r, s.prune)
  
  if (sum(r) > 0) {
    # Save selected ASVs to CSV
    pm <- psmelt(selected_phyloseq)
    write.csv(pm, file = output_csv, row.names = FALSE)
    cat("Selected ASVs saved to:", output_csv, "\n")
    return(selected_phyloseq)
  } else {
    print("No ASVs selected.")
    return(NULL)
  }
}

# Example usage:
#rf_physeq <- RandomForest_selected_ASVs(physeq,response_var = "host.species",na_vars = c("Habitat", "Ecoregion_III", "Host_genus", "Age"))

