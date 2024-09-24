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
#' @param n_top_predictors An integer specifying the number of top predictors to select. Default is 50.
#' @param output_csv A character string specifying the output CSV file name. Default is "randomforest_results.csv".
#' @param na_vars A character vector specifying the sample data variables to check for NA values. Default is NULL.
#' @return A pruned phyloseq object containing the selected ASVs.
#' @importFrom phyloseq subset_taxa sample_data tax_table sample_sums prune_samples prune_taxa otu_table psmelt taxa_sums
#' @importFrom dplyr arrange desc
#' @importFrom randomForest randomForest importance
#' @importFrom utils write.csv
#' @examples
#' if (interactive()) {
#'   # Ensure you have a phyloseq object `physeq`
#'   # Use the name of the column in sample_data for response_var
#'   rf_physeq <- RandomForest_selected_ASVs(physeq_16SOTU, response_var = "Host.genus",
#'     na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet"))
#' }
#' @export
RandomForest_selected_ASVs <- function(physeq, response_var, minlib = 15000, prunescale = 0.0001, ntree = 100, n_top_predictors = 50, output_csv = "randomforest_results.csv", na_vars = NULL) {
  suppressMessages({
    
    # Check if response_var exists in sample_data
    if (!response_var %in% colnames(phyloseq::sample_data(physeq))) {
      stop(paste("Response variable", response_var, "not found in sample data."))
    }
    
    # Ensure na_vars includes response_var
    if (!is.null(na_vars)) {
      na_vars <- unique(c(response_var, na_vars))
    } else {
      na_vars <- response_var
    }
    
    # Remove samples with NA in the specified variables
    sample_data_df <- data.frame(phyloseq::sample_data(physeq))
    for (var in na_vars) {
      sample_data_df <- sample_data_df[!is.na(sample_data_df[[var]]), ]
    }
    phyloseq::sample_data(physeq) <- phyloseq::sample_data(sample_data_df)
    
    # Convert categorical data to factors
    for (var in na_vars) {
      phyloseq::sample_data(physeq)[[var]] <- as.factor(phyloseq::sample_data(physeq)[[var]])
    }
    
    # Check if response variable has at least two classes
    if (length(unique(sample_data_df[[response_var]])) < 2) {
      stop("The response variable needs at least two classes for classification.")
    }
    
    # Remove taxa with any NA values in tax_table
    physeq <- phyloseq::subset_taxa(physeq, apply(phyloseq::tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
    
    # Remove samples with zero reads after subsetting
    physeq <- phyloseq::prune_samples(phyloseq::sample_sums(physeq) > 0, physeq)
    
    # Prune out rare OTUs by mean relative abundance set by prunescale
    tax.mean <- phyloseq::taxa_sums(physeq) / phyloseq::nsamples(physeq)
    s.prune <- phyloseq::prune_taxa(tax.mean > prunescale * minlib, physeq)
    
    # Replace empty or NA values in tax_table with "unidentified"
    phyloseq::tax_table(s.prune) <- replace(phyloseq::tax_table(s.prune), is.na(phyloseq::tax_table(s.prune)), "unidentified")
    phyloseq::tax_table(s.prune) <- replace(phyloseq::tax_table(s.prune), phyloseq::tax_table(s.prune) == "", "unidentified")
    
    predictors <- t(phyloseq::otu_table(s.prune))
    
    # Make one column for our outcome/response variable 
    response <- as.factor(phyloseq::sample_data(s.prune)[[response_var]])
    
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
    sal.classify <- randomForest::randomForest(response ~ ., data = rf.data, ntree = ntree, keep.forest = FALSE, proximity = FALSE, importance = TRUE)
    
    # Extract variable importance
    imp <- randomForest::importance(sal.classify)
    imp_df <- data.frame(predictors = rownames(imp), imp)
    imp_df <- dplyr::arrange(imp_df, dplyr::desc(MeanDecreaseGini))
    imp_df$predictors <- factor(imp_df$predictors, levels = imp_df$predictors)
    
    # Select the top n_top_predictors 
    imp_top <- imp_df[1:n_top_predictors, ]
    
    # Subset selected ASVs from phyloseq object
    otunames <- imp_top$predictors
    r <- rownames(phyloseq::tax_table(s.prune)) %in% otunames
    selected_phyloseq <- phyloseq::prune_taxa(r, s.prune)
    
    if (sum(r) > 0) {
      # Save selected ASVs to CSV
      pm <- phyloseq::psmelt(selected_phyloseq)
      utils::write.csv(pm, file = output_csv, row.names = FALSE)
      cat("Selected ASVs saved to:", output_csv, "\n")
      return(selected_phyloseq)
    } else {
      print("No ASVs selected.")
      return(NULL)
    }
  })
}

# Example usage:
# rf_physeq <- RandomForest_selected_ASVs(physeq_16SOTU,response_var = "Host.genus",
# na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet"))
# saveRDS(rf_physeq,"rf_physeq.rds")
# rf_physeq@tax_table
