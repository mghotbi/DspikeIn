#' @title Create a Regression Plot with Faceting by Range
#'
#' @description This function generates a customizable scatter plot with a linear regression line,
#' statistical equation, and facets based on a specified range variable.
#' The x and y variables are transformed using a natural log transformation (`log1p`) to handle zero values.
#'
#' @param data A data frame containing the variables to plot.
#' @param x_var A string specifying the name of the x-axis variable.
#' @param y_var A string specifying the name of the y-axis variable.
#' @param custom_range A numeric vector for defining custom ranges for the 'Percentage' column (default: c(0.1, 15, 30, 50, 75, 100)).
#' @param formula A formula for the regression equation (default: y ~ x).
#' @param plot_title A string specifying the title of the plot (default: NULL, no title will be shown if not provided).
#'
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth theme_minimal labs theme element_text scale_color_manual facet_wrap
#' @importFrom ggpubr stat_regline_equation stat_cor
#' @seealso \code{\link[ggpubr]{stat_regline_equation}}, \code{\link[ggpubr]{stat_cor}}, \code{\link[ggplot2]{facet_wrap}}
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("metadata_full", package = "DspikeIn")
#'
#'   plot_object <- regression_plot(
#'     data = metadata_full,
#'     x_var = "Observed",
#'     y_var = "Total_Reads_spiked",
#'     custom_range = c(0.1, 15, 30, 50, 75, 100)
#'   )
#'
#'   # Print the plot output
#'   print(plot_object)
#' }
#' @export
regression_plot <- function(data, x_var, y_var,
                            custom_range = c(0.1, 15, 30, 50, 75, 100),
                            formula = y ~ x, plot_title = NULL) {
  # Ensure x_var and y_var exist in the dataset
  if (!(x_var %in% colnames(data)) || !(y_var %in% colnames(data))) {
    stop("Specified x_var or y_var not found in the data.")
  }

  if (!("Percentage" %in% colnames(data))) {
    stop("Column 'Percentage' is required for creating the Range variable.")
  }

  # Ensure custom_range is valid
  if (!is.numeric(custom_range) || length(custom_range) < 2 || !all(diff(custom_range) > 0)) {
    stop("custom_range must be a numeric vector with at least two increasing values.")
  }

  # Remove NA values in relevant columns
  data <- data[!is.na(data[[x_var]]) & !is.na(data[[y_var]]) & !is.na(data$Percentage), ]

  # Apply log transformation using log1p (log(1 + x))/prev 0
  data[[x_var]] <- log1p(data[[x_var]])
  data[[y_var]] <- log1p(data[[y_var]])

  # Create the Range variable using custom ranges
  range_labels <- paste0(head(custom_range, -1), "-", tail(custom_range, -1), "%")
  data$Range <- cut(data$Percentage, breaks = custom_range, include.lowest = TRUE, labels = range_labels)

  # Remove any additional NA values from factor conversion
  data <- data[!is.na(data$Range), ]

  # Use DspikeIn color palette
  color_palette <- DspikeIn::color_palette$cool_MG
  num_ranges <- length(range_labels)

  if (length(color_palette) < num_ranges) {
    stop("Not enough colors in 'DspikeIn::color_palette$mix_MG' to match the number of range groups.")
  }

  # Generate the plot
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Range), size = 3, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    ggplot2::scale_color_manual(values = color_palette[seq_len(num_ranges)]) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("log1p(", x_var, ")"),
      y = paste0("log1p(", y_var, ")"),
      color = "Range"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold"),
      plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(size = 13, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.size = ggplot2::unit(1, "cm")
    ) +
    ggpubr::stat_regline_equation(ggplot2::aes(label = after_stat(eq.label)),
      formula = formula,
      label.x.npc = "left", label.y.npc = 0.85,
      vjust = 1, hjust = -0.01
    ) +
    ggpubr::stat_cor(ggplot2::aes(label = after_stat(paste(..rr.label.., ..p.label.., sep = "~~~"))),
      label.x.npc = "left", label.y.npc = 0.75
    ) +
    ggplot2::facet_wrap(~Range, scales = "free")

  # Add title if provided
  if (!is.null(plot_title)) {
    plot <- plot + ggplot2::labs(title = plot_title)
  }

  return(plot)
}


# Example usage;
# plot_object <- regression_plot(
# data = metadata,
# x_var = "Observed",  #  metadata is a data frame fromat
# y_var = "Spiked_Reads",
#  custom_range = c(0.1, 15, 30, 50, 75, 100),  # Define percentage ranges
#  plot_title = NULL  # No title by default
# )
# print(plot_object)
# 
