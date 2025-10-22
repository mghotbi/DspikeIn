#' @title Calculate Summary Statistics Table
#'
#' @description
#' Computes summary statistics (mean, standard deviation, standard error, and quartiles)
<<<<<<< HEAD
#' for numeric columns of a data frame. Generates a formatted flextable, and saves
#' both \code{.docx} and \code{.csv} versions of the table.
#'
#' @param data A data frame containing numeric variables.
#' @param output_path Optional. Character string specifying the DOCX output path.
#' Default is \code{"post_eval_summary.docx"}.
#'
#' @return A \code{flextable} object with formatted summary statistics.
#'
#' @importFrom flextable flextable fontsize font color bold italic save_as_docx
=======
#' for numeric columns of a data frame, and saves the resulting table as a `.csv` file.
#' This utility is designed for quick post-analysis summaries and ensures that
#' all numerical columns are summarized consistently across datasets.
#'
#' @param data A data frame containing numeric variables.
#' @param output_path Optional. Character string specifying the CSV output path.
#' Default is \code{"post_eval_summary.csv"}.
#'
#' @return A data frame containing mean, standard deviation, standard error,
#' first quartile (Q25), median, and third quartile (Q75) for each numeric column.
#'
#' @details
#' This function provides a concise statistical overview of numeric datasets.
#' It computes central tendency, dispersion, and quartile-based spread metrics
#' and writes the results to a `.csv` file suitable for reporting or downstream use.
#'
>>>>>>> bioc-rev1
#' @importFrom dplyr select_if summarise_all
#' @importFrom utils write.csv
#' @importFrom stats sd quantile median
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # --- Phyloseq example ---
#'   data("physeq", package = "DspikeIn")
#'   absolute_count <- phyloseq::otu_table(physeq)
#'
<<<<<<< HEAD
#'   tmp_docx <- file.path(tempdir(), "physeq_summary.docx")
#'   summary_table_physeq <- calculate_summary_stats_table(
#'     data = as.data.frame(absolute_count),
#'     output_path = tmp_docx
#'   )
#'   print(summary_table_physeq)
#'   if (file.exists(tmp_docx)) file.remove(tmp_docx)
=======
#'   tmp_csv <- file.path(tempdir(), "physeq_summary.csv")
#'   summary_table_physeq <- calculate_summary_stats_table(
#'     data = as.data.frame(absolute_count),
#'     output_path = tmp_csv
#'   )
#'   print(summary_table_physeq)
#'   if (file.exists(tmp_csv)) file.remove(tmp_csv)
>>>>>>> bioc-rev1
#'
#'   # --- TSE example ---
#'   data("tse", package = "DspikeIn")
#'   tse_counts <- SummarizedExperiment::assay(tse)
#'
<<<<<<< HEAD
#'   tmp_docx2 <- file.path(tempdir(), "tse_summary.docx")
#'   summary_table_tse <- calculate_summary_stats_table(
#'     data = as.data.frame(tse_counts),
#'     output_path = tmp_docx2
#'   )
#'   print(summary_table_tse)
#'   if (file.exists(tmp_docx2)) file.remove(tmp_docx2)
=======
#'   tmp_csv2 <- file.path(tempdir(), "tse_summary.csv")
#'   summary_table_tse <- calculate_summary_stats_table(
#'     data = as.data.frame(tse_counts),
#'     output_path = tmp_csv2
#'   )
#'   print(summary_table_tse)
#'   if (file.exists(tmp_csv2)) file.remove(tmp_csv2)
>>>>>>> bioc-rev1
#' }
#'
#' @export
calculate_summary_stats_table <- function(data, output_path = NULL) {
<<<<<<< HEAD
  if (!is.data.frame(data)) stop("'data' must be a data frame.")

  numeric_data <- dplyr::select_if(data, is.numeric)

  if (ncol(numeric_data) == 0) stop("No numeric columns found in the data.")

  # Calculate statistics
  summary_stats <- dplyr::summarise_all(numeric_data, list(
    mean = ~ mean(., na.rm = TRUE),
    sd = ~ sd(., na.rm = TRUE),
    se = ~ sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
    q25 = ~ quantile(., 0.25, na.rm = TRUE),
    median = ~ median(., na.rm = TRUE),
    q75 = ~ quantile(., 0.75, na.rm = TRUE)
  ))

  # Format
  ft <- flextable::flextable(summary_stats) |>
    flextable::fontsize(size = 10) |>
    flextable::font(part = "all", fontname = "Inconsolata") |>
    flextable::color(part = "header", color = "#582F0E") |>
    flextable::bold(part = "header") |>
    flextable::italic()

  if (is.null(output_path)) output_path <- "post_eval_summary.docx"

  # Save DOCX
  flextable::save_as_docx(ft, path = output_path)

  # Save CSV
  csv_path <- sub("\\.docx$", ".csv", output_path)
  utils::write.csv(summary_stats, file = csv_path, row.names = FALSE)

  message("DOCX table saved to: ", output_path)
  message("CSV summary saved to: ", csv_path)

  return(ft)
}

# Example usage:
# summary_table <- calculate_summary_stats_table(absolute_counts, output_path = "summary.docx")
=======
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  
  # Select only numeric columns
  numeric_data <- dplyr::select_if(data, is.numeric)
  if (ncol(numeric_data) == 0) {
    stop("No numeric columns found in the data.")
  }
  
  # Calculate summary statistics
  summary_stats <- dplyr::summarise_all(numeric_data, list(
    mean = ~ mean(., na.rm = TRUE),
    sd = ~ stats::sd(., na.rm = TRUE),
    se = ~ stats::sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
    q25 = ~ stats::quantile(., 0.25, na.rm = TRUE),
    median = ~ stats::median(., na.rm = TRUE),
    q75 = ~ stats::quantile(., 0.75, na.rm = TRUE)
  ))
  
  # Save as CSV
  if (is.null(output_path)) {
    output_path <- "post_eval_summary.csv"
  }
  utils::write.csv(summary_stats, file = output_path, row.names = FALSE)
  
  message("CSV summary saved to: ", normalizePath(output_path, mustWork = FALSE))
  
  return(summary_stats)
}

# Example usage:
# summary_table <- calculate_summary_stats_table(absolute_counts, output_path = "summary.csv")
>>>>>>> bioc-rev1
# print(summary_table)
# post_eval_summary <- read.csv("summary.csv")
