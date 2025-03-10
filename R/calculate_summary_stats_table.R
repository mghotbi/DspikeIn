#' @title Calculate Summary Statistics Table
#'
#' @description This function calculates summary statistics (mean, standard deviation, standard error, quantiles) for numeric columns in a data frame,
#' creates a flextable for formatted output, and saves the table in a Word document and CSV file.
#'
#' @param data A data frame containing the data to be summarized.
#' @param output_path A character string specifying the output path for the Word document. Default is NULL, which sets the output path to "post_eval_summary.docx".
#' @return A flextable object containing the summary statistics.
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE) &&
#'  requireNamespace("microbiome", quietly = TRUE)) {
#'   # Load necessary libraries
#'   library(microbiome)
#'
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Extract absolute counts using the microbiome package
#'   absolute_count <- microbiome::meta(physeq_16SOTU)
#'
#'   # Calculate summary statistics and save the table to a Word document
#'   summary_table <- calculate_summary_stats_table(absolute_count, output_path = "summary.docx")
#'
#'   # Print the calculated summary table
#'   print(summary_table)
#' }
#' }
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
#' @importFrom dplyr select_if summarise_all
#' @importFrom utils write.csv
#' @export
calculate_summary_stats_table <- function(data, output_path = NULL) {
  suppressMessages({

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
      flextable::color(part = "header", color = "#582F0E") %>%
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
    cat("\U0001F4BE	Table saved in docx format:", output_path, "\n")
    cat("\U0001F4BE	Summary statistics saved as CSV:", csv_path, "\n")

    return(ft)
  })
}

# Example usage:
# summary_table <- calculate_summary_stats_table(absolute_counts, output_path = "summary.docx")
# print(summary_table)
# post_eval_summary <- read.csv("summary.csv")
