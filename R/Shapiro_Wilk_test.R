#' Perform Shapiro-Wilk Test for Normality of Residuals with Multiple Transformations
#'
#' This function performs the Shapiro-Wilk test for normality of residuals on numeric variables
#' within a data frame, using specified transformations. It prints the p-value for each variable
#' and transformation combination.
#'
#' @param data A data frame containing the data to test.
#' @param grouping_var A character string specifying the column name to use for grouping.
#' @return NULL. The function prints the Shapiro-Wilk test results for each variable and transformation.
#' @examples
#' # Example usage:
#' # Assume `methods_df` is a data frame with a column "Methods" for grouping
#' # and other numeric columns to test.
#' # Shapiro_Wilk_test(methods_df, "Methods")
#' @export
Shapiro_Wilk_test <- function(data, grouping_var) {
  # Ensure the grouping variable is a factor
  data[[grouping_var]] <- as.factor(data[[grouping_var]])
  
  # Select numeric variables
  numeric_vars <- dplyr::select_if(data, is.numeric)
  
  for (var_name in names(numeric_vars)) {
    clean_var <- na.omit(numeric_vars[[var_name]])
    
    # Shapiro-Wilk test for normality of residuals without transformation
    shapiro_test_raw <- shapiro.test(resid(aov(clean_var ~ data[[grouping_var]])))
    print(paste("Shapiro-Wilk test for", var_name, "without transformation: p-value =", shapiro_test_raw$p.value))
    
    # Handle square root transformation (remove non-positive values)
    sqrt_var <- clean_var[clean_var >= 0]
    if (length(sqrt_var) > 1) {
      shapiro_test_sqrt <- shapiro.test(resid(aov(sqrt(sqrt_var) ~ data[[grouping_var]][clean_var >= 0])))
      print(paste("Shapiro-Wilk test for", var_name, "with square root transformation: p-value =", shapiro_test_sqrt$p.value))
    } else {
      print(paste("Not enough positive values for square root transformation of", var_name))
    }
    
    # Handle log transformation (remove non-positive values)
    log_var <- clean_var[clean_var > 0]
    if (length(log_var) > 1) {
      shapiro_test_log <- shapiro.test(resid(aov(log(log_var) ~ data[[grouping_var]][clean_var > 0])))
      print(paste("Shapiro-Wilk test for", var_name, "with log transformation: p-value =", shapiro_test_log$p.value))
    } else {
      print(paste("Not enough positive values for log transformation of", var_name))
    }
    
    # Handle z-scale transformation
    if (length(clean_var) > 1) {
      shapiro_test_zscale <- shapiro.test(resid(aov(scale(clean_var) ~ data[[grouping_var]])))
      print(paste("Shapiro-Wilk test for", var_name, "with Z-scaling: p-value =", shapiro_test_zscale$p.value))
    } else {
      print(paste("Not enough values for Z-scaling of", var_name))
    }
  }
}

# Example:
#Shapiro_Wilk_test(methods, "Methods")

