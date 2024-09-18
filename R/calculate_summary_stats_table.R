#' Calculate Summary Statistics Table
#'
#' This function calculates summary statistics (mean, standard deviation, standard error, quantiles) for numeric columns in a data frame,
#' creates a flextable for formatted output, and saves the table in a Word document and CSV file.
#'
#' @param data A data frame containing the data to be summarized.
#' @param output_path A character string specifying the output path for the Word document. Default is NULL, which sets the output path to "post_eval_summary.docx".
#' @return A flextable object containing the summary statistics.
#' @examples
#' \dontrun{
#' # Calculate summary statistics and save the table to a Word document
#' summary_table <- calculate_summary_stats_table(initial_stat_ASV, output_path = "summary.docx")
#' 
#' # Read a post-evaluation summary from a CSV file
#' post_eval_summary <- read.csv("post_eval_summary.csv")
#'
#' # Print the calculated summary table
#' print(summary_table)
#'
#' # Preview the post-evaluation summary
#' head(post_eval_summary)
#' }
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @export
calculate_summary_stats_table <- function(data, output_path = NULL) {
  suppressMessages({
    # Load necessary libraries
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required but not installed.")
    }
    if (!requireNamespace("flextable", quietly = TRUE)) {
      stop("Package 'flextable' is required but not installed.")
    }
    
    # Calculate summary statistics for numeric columns
    summary_stats <- dplyr::select_if(data, is.numeric) %>%
      dplyr::summarise_all(list(
        mean = ~mean(.),
        sd = ~ifelse(all(is.na(.)), NA, stats::sd(., na.rm = TRUE)),
        se = ~ifelse(all(is.na(.)), NA, stats::sd(., na.rm = TRUE) / sqrt(length(.))),
        q25 = ~stats::quantile(., 0.25, na.rm = TRUE),
        median = ~stats::median(.),
        q75 = ~stats::quantile(., 0.75, na.rm = TRUE)
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
    utils::write.csv(summary_stats, file = csv_path, row.names = FALSE)
    
    # Print a message indicating where the files were saved
    cat("Table saved in docx format:", output_path, "\n")
    cat("Summary statistics saved as CSV:", csv_path, "\n")
    
    return(ft)
  })
}

# Example usage:
# post_eval_summary <- calculate_summary_stats_table(absolute_counts)
# print(post_eval_summary)
# summary_table <- calculate_summary_stats_table(initial_stat_ASV, output_path = "summary.docx")
# post_eval_summary <- read.csv("summary.csv")
