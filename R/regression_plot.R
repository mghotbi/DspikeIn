#' Create a Regression Plot with Faceting by Range
#'
#' This function generates a customizable scatter plot with a linear regression line, 
#' statistical equation, and facets based on a specified range variable.
#'
#' @param data A data frame containing the variables to plot.
#' @param x_var A string specifying the name of the x-axis variable.
#' @param y_var A string specifying the name of the y-axis variable.
#' @param custom_range A numeric vector for defining custom ranges for the 'Percentage' column (default: c(0.1, 15, 30, 50, 75, 100)).
#' @param formula A formula for the regression equation (default: y ~ x).
#' @param plot_title A string specifying the title of the plot (default: NULL, no title will be shown if not provided).
#' @param pseudocount A small numeric value that is added to both the x and y variables
#' to avoid issues when there are zero values in the data. This is helpful in situations 
#' where zero values would cause issues with calculations, such as log transformations 
#' or statistical modeling. Default is 1e-6.
#'
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth theme_minimal labs theme element_text scale_color_manual facet_wrap
#' @importFrom ggpubr stat_regline_equation stat_cor
#' 
#' @examples
#' \dontrun{
#' # Example of how to call the function with custom ranges
#' plot_object <- regression_plot(
#'   data = metadata, 
#'   x_var = "Richness", 
#'   y_var = "Total_Reads_spiked", 
#'   custom_range = c(0.1, 15, 30, 50, 75, 100)
#' )
#' print(plot_object)
#' }
#' @export
regression_plot <- function(data, x_var, y_var, custom_range = c(0.1, 15, 30, 50, 75, 100), 
                            formula = y ~ x, plot_title = NULL, pseudocount = 1e-6) {
  
  # Ensure x_var and y_var exist in the dataset
  if (!(x_var %in% colnames(data)) || !(y_var %in% colnames(data))) {
    stop("Specified x_var or y_var not found in the data.")
  }
  
  # Add a pseudocount to avoid issues with zero values
  data[[x_var]] <- data[[x_var]] + pseudocount
  data[[y_var]] <- data[[y_var]] + pseudocount
  
  # Check if Percentage column exists for faceting
  if (!("Percentage" %in% colnames(data))) {
    stop("Column 'Percentage' is required for creating the Range variable.")
  }
  
  # Create the Range variable using custom ranges
  data$Range <- cut(data$Percentage, breaks = custom_range, include.lowest = TRUE, 
                          labels = paste0(custom_range[-length(custom_range)], "-", custom_range[-1], "%"))
  
  # Basic ggplot
  plot <- ggplot2::ggplot(data, ggplot2::aes_string(x = x_var, y = y_var)) +
    ggplot2::geom_point(ggplot2::aes_string(color = "Range"), size = 3, alpha = 0.7) +  # Scatter plot with color based on Range
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +  # Regression line
    ggplot2::scale_color_manual(values = c(  # Custom colors for Range
      "0.1-15%" = "#FF7F00", "15-30%" = "#E41A1C", "30-50%" = "firebrick4", 
      "50-75%" = "#4daf4a", "75-100%" = "#984EA3"
    )) +
    ggplot2::theme_minimal() +  # Clean theme
    ggplot2::labs(
      x = x_var,  # X-axis label
      y = y_var,  # Y-axis label
      color = "Range"  # Legend title
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
    ggpubr::stat_regline_equation(ggplot2::aes(label = ..eq.label..), formula = formula, 
                                  label.x.npc = "left", label.y.npc = 0.85, 
                                  vjust = 1, hjust = -0.01) +  # Regression equation closer
    ggpubr::stat_cor(ggplot2::aes(label = paste(..rr.label.., ..p.label.., sep = "~~~")), 
                     label.x.npc = "left", label.y.npc = 0.75) +  
    ggplot2::facet_wrap(~ Range, scales = "free")  
  
  # Add title only if plot_title is not NULL
  if (!is.null(plot_title)) {
    plot <- plot + ggplot2::labs(title = plot_title)
  }
  
  return(plot)
}

# Example usage;
#plot_object <- regression_plot(
#data = metadata, 
# x_var = "Richness.x",  #  metadata is a data frame fromat
# y_var = "Total_Reads_spiked",  
#  custom_range = c(0.1, 15, 30, 50, 75, 100),  # Define percentage ranges
#  plot_title = NULL  # No title by default
#)
# print(plot_object)
