#' @title Perform Shapiro-Wilk Test for Normality with Transformations
#' @description This function performs the Shapiro-Wilk test for normality of residuals
#' on numeric variables within a data frame, using different transformations.
#'
#' @param data A data frame containing the data to test.
#' @param grouping_var A character string specifying the column name to use for grouping.
#' @return NULL. The function prints the Shapiro-Wilk test results for each variable and transformation.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("methods", package = "DspikeIn")
#'   # Perform the Shapiro-Wilk test on the dataset
#'   Shapiro_Wilk_test(methods, "Methods")
#' }
#' }
#' @importFrom stats shapiro.test aov resid
#' @importFrom dplyr select_if
#' @export
Shapiro_Wilk_test <- function(data, grouping_var) {
  suppressMessages({
    # Validate input
    if (!is.data.frame(data)) stop("\U0000274C Error: `data` must be a data frame.")
    if (!grouping_var %in% colnames(data)) stop("\U0000274C Error: `grouping_var` not found in data.")

    # Ensure the grouping variable is a factor
    data[[grouping_var]] <- as.factor(data[[grouping_var]])

    # Select numeric variables
    numeric_vars <- dplyr::select_if(data, is.numeric)

    for (var_name in names(numeric_vars)) {
      clean_var <- na.omit(numeric_vars[[var_name]])

      message("\n **Testing Normality for**: ", var_name)

      # Perform Shapiro-Wilk test (raw values)
      shapiro_test_raw <- stats::shapiro.test(stats::resid(stats::aov(clean_var ~ data[[grouping_var]])))
      message(" \U00002139  No Transformation: p-value = ", formatC(shapiro_test_raw$p.value, format = "e", digits = 3))

      # Square root transformation (only for non-negative values)
      sqrt_var <- clean_var[clean_var >= 0]
      if (length(sqrt_var) > 1) {
        shapiro_test_sqrt <- stats::shapiro.test(stats::resid(stats::aov(sqrt(sqrt_var) ~ data[[grouping_var]][clean_var >= 0])))
        message(" \U00002139  Square Root Transformation: p-value = ", formatC(shapiro_test_sqrt$p.value, format = "e", digits = 3))
      } else {
        message("\U00002139 Not enough positive values for square root transformation.")
      }

      # Log transformation (only for positive values)
      log_var <- clean_var[clean_var > 0]
      if (length(log_var) > 1) {
        shapiro_test_log <- stats::shapiro.test(stats::resid(stats::aov(log(log_var) ~ data[[grouping_var]][clean_var > 0])))
        message(" \U00002139  Log Transformation: p-value = ", formatC(shapiro_test_log$p.value, format = "e", digits = 3))
      } else {
        message(" \U00002139 Not enough positive values for log transformation.")
      }

      # Z-score standardization
      if (length(clean_var) > 1) {
        shapiro_test_zscale <- stats::shapiro.test(stats::resid(stats::aov(scale(clean_var) ~ data[[grouping_var]])))
        message(" \U00002139 Z-Scaling: p-value = ", formatC(shapiro_test_zscale$p.value, format = "e", digits = 3))
      } else {
        message(" \U00002139  Not enough values for Z-scaling.")
      }
    }
  })
}

# Example Usage:
# file_path <- system.file("extdata", "methods.rds", package = "DspikeIn")
# methods <- readRDS(file_path)
# Shapiro_Wilk_test(methods, "Methods")
