#' Calculate Summary Statistics Table
#'
#' This function calculates summary statistics (mean, standard deviation, standard error, quantiles) for numeric columns in a data frame,
#' creates a flextable for formatted output, and saves the table in a Word document and CSV file.
#'
#' @param data A data frame containing the data to be summarized.
#' @param output_path A character string specifying the output path for the Word document. Default is NULL, which sets the output path to "post_eval_summary.docx".
#' @return A flextable object containing the summary statistics.
#' @examples
#' # Calculate summary statistics and save the table
#' summary_table <- calculate_summary_stats_table(initial_stat_ASV, output_path = "summary.docx")
#' post_eval_summary <- read.csv("post_eval_summary.csv")
#' @export
calculate_summary_stats_table <- function(data, output_path = NULL) {
  library(dplyr)
  library(flextable)
  
  # Calculate summary statistics for numeric columns
  summary_stats <- data %>%
    group_by() %>%
    select_if(is.numeric) %>%
    summarise_all(list(
      mean = ~mean(.),
      sd = ~ifelse(all(is.na(.)), NA, sd(., na.rm = TRUE)),
      se = ~ifelse(all(is.na(.)), NA, sd(., na.rm = TRUE) / sqrt(length(.))),
      q25 = ~quantile(., 0.25, na.rm = TRUE),
      median = ~median(.),
      q75 = ~quantile(., 0.75, na.rm = TRUE)
    ))
  
  # Create a flextable for the summary statistics
  ft <- flextable::flextable(summary_stats) %>% 
    flextable::fontsize(size = 10) %>% 
    flextable::font(part = "all", fontname = "Inconsolata") %>% 
    flextable::color(part = "header", color = "red4") %>% 
    flextable::bold(part = "header") %>% 
    flextable::italic()
  
  # Set default output directory if none provided
  if (is.null(output_path)) {
    output_path <- "post_eval_summary.docx"
  }
  
  # Save the flextable as a Word document
  flextable::save_as_docx(ft, path = output_path)
  
  # Save summary statistics data frame as CSV
  csv_path <- sub(".docx", ".csv", output_path)
  write.csv(summary_stats, file = csv_path, row.names = FALSE)
  
  # Print a message indicating where the files were saved
  cat("Table saved in docx format:", output_path, "\n")
  cat("Summary statistics saved as CSV:", csv_path, "\n")
  
  return(ft)
}

# Example usage:
# Calculate summary statistics and save the table
# summary_table <- calculate_summary_stats_table(initial_stat_ASV, output_path = "summary.docx")
# post_eval_summary <- read.csv("summary.csv")
