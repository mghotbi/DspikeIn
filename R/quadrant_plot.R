#' @title Generate Custom Quadrant Plots for Node Metrics
#'
#' @description This function visualizes relationships between **any two node metrics** in a quadrant plot.
#' Quadrant labels dynamically adjust based on the selected X and Y axis metrics.
#'
#' @param metrics A `data.frame` containing computed node metrics.
#' @param x_metric Character. Column name to use for the x-axis. Default is `"Degree"`.
#' @param y_metric Character. Column name to use for the y-axis. Default is `"Redundancy"`.
#' @param x_threshold Numeric. X-axis threshold for quadrant separation. Default is `median(x_metric)`.
#' @param y_threshold Numeric. Y-axis threshold for quadrant separation. Default is `median(y_metric)`.
#' @param top_quantile Numeric. Quantile threshold (0-1) to highlight top nodes. Default is `0.95`.
#' @param point_size Numeric. Size of the points. Default is `3`.
#'
#' @return A `ggplot` object representing the customized quadrant plot.
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   g <- load_graphml("Complete.graphml")
#'
#'   # Compute node-level metrics
#'   result <- node_level_metrics(g)
#'   metrics <- result$metrics
#'
#'   # Generate a quadrant plot using Degree and Efficiency
#'   plot <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Efficiency")
#'   print(plot)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_color_manual theme_minimal labs
#' @importFrom dplyr filter mutate case_when
#' @importFrom ggrepel geom_text_repel
#' @export
quadrant_plot <- function(metrics,
                          x_metric = "Degree",
                          y_metric = "Redundancy",
                          x_threshold = NULL,
                          y_threshold = NULL,
                          top_quantile = 0.95,
                          point_size = 3) {
  message("\U0001F578 Generating customizable network visualization...")

  # Validate input
  if (!(x_metric %in% names(metrics)) || !(y_metric %in% names(metrics))) {
    stop("\U0000274C Error: Specified x_metric or y_metric not found in metrics data frame.")
  }

  # Compute Default Thresholds If Not Provided
  if (is.null(x_threshold)) x_threshold <- median(metrics[[x_metric]], na.rm = TRUE)
  if (is.null(y_threshold)) y_threshold <- median(metrics[[y_metric]], na.rm = TRUE)

  # Assign Quadrant Labels Dynamically Based on Metrics
  metrics <- metrics %>%
    mutate(
      Quadrant = case_when(
        .data[[x_metric]] >= x_threshold & .data[[y_metric]] >= y_threshold ~ paste("High", y_metric, "& High", x_metric),
        .data[[x_metric]] >= x_threshold & .data[[y_metric]] < y_threshold ~ paste("Low", y_metric, "& High", x_metric),
        .data[[x_metric]] < x_threshold & .data[[y_metric]] >= y_threshold ~ paste("High", y_metric, "& Low", x_metric),
        TRUE ~ paste("Low", y_metric, "& Low", x_metric)
      )
    )

  # Define Custom Colors for Quadrants (Dynamic)
  quadrant_colors <- setNames(
    c("#50A5B1", "#F8D49B", "#9C8DCB", "#E06C75"),
    unique(metrics$Quadrant)
  )

  # Identify Key Nodes to Highlight (Top Quantile)
  top_nodes <- dplyr::filter(metrics, .data[[y_metric]] > quantile(metrics[[y_metric]], top_quantile, na.rm = TRUE))

  #  **Quadrant Plot: Custom Metric Relationships**
  plot <- ggplot2::ggplot(metrics, ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]], color = Quadrant)) +
    ggplot2::geom_point(alpha = 0.9, size = point_size) +
    ggplot2::scale_color_manual(values = quadrant_colors) +
    ggrepel::geom_text_repel(data = top_nodes, ggplot2::aes(label = Node), size = 3, color = "black", max.overlaps = 25) +
    ggplot2::geom_vline(xintercept = x_threshold, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = y_threshold, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = paste("Quadrant Analysis:", x_metric, "vs.", y_metric),
      x = x_metric,
      y = y_metric,
      color = "Quadrant"
    ) +
    ggplot2::theme_minimal(base_size = 13)

  return(plot)
}


# Usage Example:
# Compute Node-Level Metrics
# result <- node_level_metrics(Complete)
# metrics <- result$metrics
# Generate Custom Plots
# plot1 <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Redundancy")
# plot2 <- quadrant_plot(metrics, x_metric = "Degree", y_metric = "Efficiency")
# plot3 <- quadrant_plot(metrics, x_metric = "Within_Module_Connectivity", y_metric = "Betweenness")
# plot4 <- quadrant_plot(metrics, x_metric = "Among_Module_Connectivity", y_metric = "Degree")

# print(plot1)
# print(plot2)
# print(plot3)
# print(plot4)
