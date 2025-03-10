#' @title Analyze and Visualize a Microbial Network
#'
#' @description Loads a microbial network from the **DspikeIn** package (`Complete.graphml`, `NoHubs.graphml`, or `NoBasid.graphml`)
#' or from a user-provided **GraphML** file. Computes network metrics, assigns modularity-based node colors (`cool_MG` from **DspikeIn**),
#' and visualizes the network using `ggraph`. Optionally saves computed **global network metrics** as a CSV file.
#'
#' @param graph_path Character. The **GraphML file path**. Default (`NULL`) loads `"Complete.graphml"` from **DspikeIn**.
#' Valid internal options: `"Complete.graphml"`, `"NoHubs.graphml"`, `"NoBasid.graphml"`.
#' If an **external path** is provided, the function loads that instead.
#'
#' @param save_metrics Logical. If `TRUE`, saves `"Global_Network_Metrics.csv"`.
#'
#' @return A **list** containing:
#'   \item{plot}{A **ggplot2** object displaying the network.}
#'   \item{metrics}{A **data.frame** with global network metrics.}
#'   \item{graph}{The **igraph** object with assigned attributes.}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   Complete <- load_graphml("Complete.graphml")
#'
#'   # Run network weighting function on the loaded dataset
#'   result <- weight_Network()
#'
#'   # Load a specific GraphML dataset from DspikeIn
#'   result <- weight_Network(graph_path = "NoHubs.graphml")
#'
#'   # Load a custom GraphML file from user directory,
#'   # for external graphml please use **full address**
#'   print(result$plot)
#'
#'   # View network metrics
#'   result$metrics
#' }
#' }
#'
#' @importFrom igraph cluster_fast_greedy membership degree edge_density is_connected
#' @importFrom igraph components induced_subgraph V E modularity transitivity diameter mean_distance vcount ecount
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point scale_edge_colour_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales rescale
#' @importFrom ggplot2 theme_void labs theme scale_color_manual element_rect guides guide_legend aes
#' @importFrom utils write.csv
#' @export
weight_Network <- function(graph_path = NULL, save_metrics = TRUE) {

  # =====================
  # Load GraphML File
  # =====================
  default_graphs <- c("Complete.graphml", "NoHubs.graphml", "NoBasid.graphml")

  if (is.null(graph_path) || graph_path %in% default_graphs) {
    # Load internal GraphML file from package
    if (is.null(graph_path)) graph_path <- "Complete.graphml"
    message("\U0001F4C2 Loading default GraphML file: ", graph_path)
    graph <- load_graphml(graph_path)
  } else {
    # Load user-provided GraphML file
    if (!file.exists(graph_path)) stop("\U0000274C Error: GraphML file not found: ", graph_path)
    message("\U0001F4C1 Loading external GraphML file: ", graph_path)
    graph <- load_graphml(graph_path)
  }

  # =====================
  # Compute Node Degree (For Node Size)
  # =====================
  igraph::V(graph)$degree <- igraph::degree(graph)

  # =====================
  # Preserve & Modify Edge Weights
  # =====================
  if (!is.null(igraph::E(graph)$weight)) {
    igraph::E(graph)$original_weight <- igraph::E(graph)$weight  # Store original weight for visualization
    igraph::E(graph)$weight <- abs(igraph::E(graph)$weight)  # Convert all weights to positive for layout calculations
  } else {
    igraph::E(graph)$weight <- rep(1, igraph::ecount(graph))
    igraph::E(graph)$original_weight <- rep(1, igraph::ecount(graph))  # Keep consistent for visualization
  }

  # =====================
  # Compute Modularity-Based Colors
  # =====================
  modules <- as.factor(igraph::membership(igraph::cluster_fast_greedy(graph)))
  node_colors <- rep(DspikeIn::color_palette$cool_MG, length.out = length(unique(modules)))

  # =====================
  # Generate Network Layout & Assign Coordinates
  # =====================
  message("\U0001F5A7 Generating network layout...")
  layout_data <- ggraph::create_layout(graph, layout = "fr")

  # Ensure x and y coordinates are assigned safely
  layout_x <- layout_data$x
  layout_y <- layout_data$y

  # =====================
  # Compute Global Metrics
  # =====================
  message("\U0001F578 Computing global network metrics...")

  global_metrics <- data.frame(
    Nodes = igraph::vcount(graph),
    Edges = igraph::ecount(graph),
    Modularity = tryCatch(igraph::modularity(igraph::cluster_fast_greedy(graph)), error = function(e) NA),
    Density = tryCatch(igraph::edge_density(graph), error = function(e) NA),
    Transitivity = tryCatch(igraph::transitivity(graph), error = function(e) NA),
    Diameter = tryCatch(if (igraph::is_connected(graph)) igraph::diameter(graph) else NA, error = function(e) NA),
    Avg_Path_Length = tryCatch(if (igraph::is_connected(graph)) igraph::mean_distance(graph) else NA, error = function(e) NA),
    Avg_Degree = tryCatch(mean(igraph::degree(graph)), error = function(e) NA)
  )

  if (save_metrics) {
    utils::write.csv(global_metrics, "Global_Network_Metrics.csv", row.names = FALSE)
    message("\U0001F5C2 Global network metrics saved as 'Global_Network_Metrics.csv'.")
  }

  # =====================
  # Generate Network Plot
  # =====================
  message("\U0001F5A7 Generating network plot...")

  edge_sign <- factor(ifelse(igraph::E(graph)$original_weight > 0, "Positive", "Negative"))
  edge_colors <- c("Positive" = "#669BBC", "Negative" =   "#EE9B00")
  edge_thickness <- scales::rescale(abs(igraph::E(graph)$original_weight), to = c(0.3, 5))  # Use original weight to reflect true weight

  plot <- ggraph::ggraph(layout_data) +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = edge_thickness, color = edge_sign), alpha = 0.2) +
    ggraph::scale_edge_colour_manual(values = edge_colors) +
    ggraph::geom_node_point(ggplot2::aes(x = layout_x, y = layout_y, size = degree, color = modules)) +
    ggrepel::geom_text_repel(ggplot2::aes(x = layout_x, y = layout_y, label = igraph::V(graph)$name), size = 3, max.overlaps = 10) +
    ggplot2::scale_color_manual(values = node_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",    # Move legend to bottom
      legend.direction = "horizontal", # Make legend horizontal
      legend.box = "horizontal",     # Align legend elements horizontally
      legend.text = ggplot2::element_text(size = 10, face = "plain"),  # Make text readable
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.key.size = ggplot2::unit(0.5, "cm"), # Reduce legend item size
      legend.spacing.x = ggplot2::unit(0.3, "pt") # Reduce horizontal space between legend items
    ) +
      ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, keywidth = 0.5, keyheight = 0.5)
    )


  return(list(plot = plot, metrics = global_metrics, graph = graph))
}

# ================
#  Run Function with Example Files
# ================

# Load the default DspikeIn dataset (Complete.graphml)
# Complete <- load_graphml("Complete.graphml")

#result <- weight_Network()

# Load a different built-in dataset (NoHubs.graphml)
#result <- weight_Network(graph_path = "NoBasid.graphml")

# Load a user-provided GraphML file
#result <- weight_Network(graph_path = "~/network.graphml")

# View the network plot
#print(result$plot)

# Check global network metrics
#result$metrics



