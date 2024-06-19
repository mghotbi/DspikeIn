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
#' # Assume `methods_df` is a data frame with a column "Methods" for grouping
#' # and other numeric columns to test.
#' # Bartlett_test(methods_df, "Methods")
#' @export
Bartlett_test <- function(data, grouping_var, transformations = c("raw", "sqrt", "log", "zscale")) {
  # Load necessary library
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  
  library(dplyr)
  
  with(data, {
    numeric_vars <- dplyr::select_if(data, is.numeric)
    for (var_name in names(numeric_vars)) {
      clean_var <- na.omit(numeric_vars[[var_name]])
      grouping <- data[[grouping_var]][!is.na(numeric_vars[[var_name]])]
      
      for (transformation in transformations) {
        # Apply specified transformation
        if (transformation == "sqrt") {
          transformed_var <- clean_var[clean_var >= 0]
          grouping_transformed <- grouping[clean_var >= 0]
          if (length(transformed_var) > 1) {
            transformed_var <- sqrt(transformed_var)
          } else {
            transformed_var <- NA
          }
        } else if (transformation == "log") {
          transformed_var <- clean_var[clean_var > 0]
          grouping_transformed <- grouping[clean_var > 0]
          if (length(transformed_var) > 1) {
            transformed_var <- log(transformed_var)
          } else {
            transformed_var <- NA
          }
        } else if (transformation == "zscale") {
          if (length(clean_var) > 1) {
            transformed_var <- scale(clean_var)
            grouping_transformed <- grouping
          } else {
            transformed_var <- NA
          }
        } else {
          transformed_var <- clean_var  # No transformation/raw
          grouping_transformed <- grouping
        }
        
        if (!any(is.na(transformed_var))) {
          # Bartlett test for homogeneity of variances
          bartlett_result <- bartlett.test(transformed_var ~ as.factor(grouping_transformed))
          
          print(paste("Bartlett test for", var_name, "with", transformation, "transformation: p-value =", bartlett_result$p.value))
        } else {
          print(paste("Not enough valid values for", var_name, "with", transformation, "transformation"))
        }
      }
    }
  })
}

# Example:
# Bartlett_test(methods_df, "Methods")
