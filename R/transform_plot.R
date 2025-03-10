#' @title Generate and Save Boxplots with Statistical Tests
#'
#' @description This function generates boxplots for multiple y variables and performs statistical
#' tests (Kruskal-Wallis or ANOVA) to compare groups. It also performs pairwise comparisons when
#' needed and saves the plots in PNG and PDF formats.
#'
#' @param data A data frame containing the data to plot.
#' @param x_var A character string specifying the column name for the x variable (categorical/factor).
#' @param y_vars A character vector specifying the column names for the y variables.
#' @param methods_var A character string specifying the column name for the grouping variable.
#' @param color_palette A character string for the color palette from `DspikeIn::color_palette$cool_MG`.
#' @param output_prefix A character string specifying the prefix for the output file names. Default is `"plot"`.
#' @param width A numeric value specifying the width of the output plot. Default is `15`.
#' @param height A numeric value specifying the height of the output plot. Default is `13`.
#' @param stat_test A character string specifying the statistical test to use (`"kruskal.test"` or `"anova"`). Default is `"kruskal.test"`.
#' @return A list containing the ggplot2 boxplot objects and the statistical comparison results.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("methods", package = "DspikeIn")
#'
#'   # Define variables for the y-axis
#'   y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")
#'
#'   # Generate and save boxplots with Kruskal-Wallis test
#'   plot_output <- transform_plot(
#'     data = methods,
#'     x_var = "Methods",
#'     y_vars = y_vars,
#'     methods_var = "Methods",  # Ensure this parameter is needed
#'     color_palette = DspikeIn::color_palette$cool_MG,
#'     stat_test = "kruskal.test"
#'   )
#'
#'   # Print the output plot
#'   print(plot_output)
#' }
#' }
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_manual labs theme_minimal element_text ggsave
#' @importFrom ggpubr stat_compare_means
#' @importFrom dplyr mutate filter across
#' @importFrom stats kruskal.test wilcox.test
#' @export
#' @export
transform_plot <- function(
    data,
    x_var,
    y_vars,
    methods_var,
    color_palette = DspikeIn::color_palette$cool_MG,
    output_prefix = "plot",
    width = 15,
    height = 13,
    stat_test = "kruskal.test"
) {

  # Ensure x_var is a factor
  data[[x_var]] <- as.factor(data[[x_var]])

  # Ensure numeric conversion for y_vars (avoid coercion errors)
  data <- data %>%
    mutate(across(all_of(y_vars), ~ suppressWarnings(as.numeric(gsub("[^0-9.]", "", .)))))

  # Remove rows with NA values
  data <- data %>%
    filter(complete.cases(across(all_of(y_vars))))

  # Define test suffix
  test_suffix <- ifelse(stat_test == "kruskal.test", "Kruskal", "ANOVA")

  # Function to create individual boxplots
  create_boxplot <- function(y_var) {
    ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
      geom_boxplot() +
      scale_fill_manual(values = color_palette) +
      labs(x = x_var, y = y_var, title = paste(y_var, "by", x_var)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 14, angle = 30, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")
      ) +
      ggpubr::stat_compare_means(method = stat_test, label = "p.signif", size = 6)
  }

  # Create and save plots
  plots <- list()

  for (y_var in y_vars) {
    plot <- create_boxplot(y_var)
    print(plot)

    # File names
    png_filename <- paste0(output_prefix, "_", y_var, "_", test_suffix, ".png")
    pdf_filename <- paste0(output_prefix, "_", y_var, "_", test_suffix, ".pdf")

    # Save plots
    ggsave(pdf_filename, plot = plot, width = width, height = height, units = "in")
    ggsave(png_filename, plot = plot, width = width, height = height, units = "in", dpi = 500)

    cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
    plots[[y_var]] <- plot
  }

  return(plots)
}

#  **Example Usage**
# y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")
# transform_plot(data = methods, x_var = "Methods", y_vars = y_vars, methods_var = "Methods",
# color_palette = DspikeIn::color_palette$cool_MG, stat_test = "kruskal.test")
