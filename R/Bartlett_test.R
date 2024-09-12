#' Perform Bartlett's Test for Homogeneity of Variances with Multiple Transformations
#'
#' This function performs Bartlett's test for homogeneity of variances on numeric variables
#' within a data frame, using specified transformations. It prints the p-value for each variable
#' and transformation combination.
#'
#' @param data A data frame containing the data to test.
#' @param grouping_var A character string specifying the column name to use for grouping.
#' @param transformations A character vector specifying the transformations to apply. Default is c("raw", "sqrt", "log", "zscale").
#' @return NULL. The function prints the Bartlett test results for each variable and transformation.
#' @examples
#' # Example usage:
#' # Assume `methods` is a data frame with a column "Methods" for grouping
#' # and other numeric columns to test.
#' # file_path <- system.file("extdata", "methods.rds", package = "DspikeIn") 
#' methods <- readRDS(file_path) 
#' print(methods) 
#' Bartlett_test(methods, "Methods")
#' @export
Bartlett_test <- function(data, grouping_var, transformations = c("raw", "sqrt", "log", "zscale")) {
  
  # Ensure 'dplyr' is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  
  # Select numeric columns from the data frame
  numeric_vars <- dplyr::select_if(data, is.numeric)
  
  # Check if numeric variables exist
  if (ncol(numeric_vars) == 0) {
    stop("No numeric variables found in the data.")
  }
  
  # Iterate over each numeric variable
  for (var_name in names(numeric_vars)) {
    clean_var <- stats::na.omit(numeric_vars[[var_name]])  # Remove NA values
    grouping <- data[[grouping_var]][!is.na(numeric_vars[[var_name]])]
    
    if (length(unique(grouping)) < 2) {
      message(paste("Grouping variable", grouping_var, "has less than two levels. Skipping variable", var_name))
      next
    }
    
    # Apply each transformation
    for (transformation in transformations) {
      
      if (transformation == "sqrt") {
        # Square root transformation
        transformed_var <- clean_var[clean_var >= 0]  # Only non-negative values allowed
        grouping_transformed <- grouping[clean_var >= 0]
        transformed_var <- if (length(transformed_var) > 1) sqrt(transformed_var) else NA
        
      } else if (transformation == "log") {
        # Log transformation
        transformed_var <- clean_var[clean_var > 0]  # Only positive values allowed
        grouping_transformed <- grouping[clean_var > 0]
        transformed_var <- if (length(transformed_var) > 1) log(transformed_var) else NA
        
      } else if (transformation == "zscale") {
        # Z-scale transformation
        if (length(clean_var) > 1) {
          transformed_var <- scale(clean_var)
          grouping_transformed <- grouping
        } else {
          transformed_var <- NA
        }
        
      } else {
        # Raw (no transformation)
        transformed_var <- clean_var
        grouping_transformed <- grouping
      }
      
      # Ensure there are enough valid values after transformation
      if (!any(is.na(transformed_var)) && length(transformed_var) > 1) {
        # Perform Bartlett's test
        bartlett_result <- stats::bartlett.test(transformed_var ~ as.factor(grouping_transformed))
        message(paste("Bartlett test for", var_name, "with", transformation, "transformation: p-value =", bartlett_result$p.value))
      } else {
        message(paste("Not enough valid values for", var_name, "with", transformation, "transformation."))
      }
    }
  }
}

# Example usage:
# file_path <- system.file("extdata", "methods.rds", package = "DspikeIn")
# methods <- readRDS(file_path)
# Bartlett_test(methods, "Methods")
