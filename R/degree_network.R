#' @title Analyze and Visualize a Microbial Network
#'
#' @description This function loads a microbial network from a GraphML file or an `igraph` object,
#' computes node degree and modularity, assigns node sizes based on degree, and visualizes the network.
#' Edge thickness is determined by weight. It also computes and optionally saves global network metrics.
#'
#' @param graph_path Character or igraph object. Path to the GraphML file or an already loaded `igraph` object.
#' @param save_metrics Logical. If TRUE, saves the global metrics as a CSV file.
#' @param metrics_path Optional. Path to save the metrics file. Default: "Global_Network_Metrics.csv".
#' @param layout_type Character. Layout algorithm: "stress" (default), "graphopt", "fr", "mds", or "kk".
#'
#' @return A list containing:
#'   \item{plot}{A ggplot object displaying the network.}
#'   \item{metrics}{A dataframe with global network metrics.}
#'   \item{layout_data}{The layout used for plotting.}
#'   \item{graph}{The annotated igraph object.}
#'
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   Complete <- load_graphml(system.file("extdata", "Complete.graphml", package = "DspikeIn"))
#'
#'   # Save metrics to a temporary file
#'   temp_metrics <- file.path(tempdir(), "Global_Network_Metrics.csv")
#'
#'   result <- degree_network(
#'     graph_path = Complete,
#'     save_metrics = TRUE,
#'     metrics_path = temp_metrics
#'   )
#'
#'   print(result$metrics)
#'   print(result$plot)
#'
#'   # Clean up temporary file
#'   unlink(temp_metrics)
#' }
#'
#' @importFrom igraph read_graph cluster_louvain membership degree edge_density is_connected components induced_subgraph V E modularity
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales rescale
#' @importFrom ggplot2 theme_void labs theme scale_color_manual element_rect unit
#' @importFrom utils write.csv
#' @export
degree_network <- function(graph_path, save_metrics = TRUE, metrics_path = "Global_Network_Metrics.csv", layout_type = "stress") {
  # Load network
  if (inherits(graph_path, "igraph")) {
    message("Using preloaded igraph object...")
    graph <- graph_path
  } else {
    if (!file.exists(graph_path)) stop("Error: Graph file not found.")
    message("Loading graph from file...")
    graph <- igraph::read_graph(graph_path, format = "graphml")
  }

  # Ensure edge weights
  message("Processing edge weights...")
  if (!is.null(igraph::E(graph)$weight)) {
    igraph::E(graph)$original_weight <- igraph::E(graph)$weight
    igraph::E(graph)$weight <- ifelse(igraph::E(graph)$weight <= 0, 1e-3, igraph::E(graph)$weight)
  } else {
    igraph::E(graph)$weight <- rep(1, igraph::ecount(graph))
  }
  edge_weights <- igraph::E(graph)$weight

  # Compute node degree
  igraph::V(graph)$degree <- igraph::degree(graph)

  # Compute communities
  community_detection <- tryCatch(igraph::cluster_louvain(graph), error = function(e) NULL)
  communities <- if (!is.null(community_detection)) {
    igraph::membership(community_detection)
  } else {
    rep(1, igraph::vcount(graph))
  }

  # Compute global metrics
  message("Computing global network metrics...")
  is_conn <- igraph::is_connected(graph)
  largest_comp <- if (!is_conn) {
    igraph::induced_subgraph(graph, which(igraph::components(graph)$membership == which.max(igraph::components(graph)$csize)))
  } else {
    graph
  }

  global_metrics <- data.frame(
    Nodes = igraph::vcount(graph),
    Edges = igraph::ecount(graph),
    Modularity = tryCatch(igraph::modularity(community_detection), error = function(e) NA),
    Density = tryCatch(igraph::edge_density(graph), error = function(e) NA),
    Transitivity = tryCatch(igraph::transitivity(graph), error = function(e) NA),
    Diameter = tryCatch(if (is_conn) igraph::diameter(largest_comp, weights = edge_weights) else NA, error = function(e) NA),
    Avg_Path_Length = tryCatch(if (is_conn) igraph::mean_distance(largest_comp, weights = edge_weights) else NA, error = function(e) NA),
    Avg_Degree = tryCatch(mean(igraph::degree(graph)), error = function(e) NA)
  )

  # Save global metrics
  if (save_metrics) {
    utils::write.csv(global_metrics, metrics_path, row.names = FALSE)
    message("Global network metrics saved as '", metrics_path, "'.")
  }

  # Safe layout generation
  safely_layout <- function(graph, layout_type, weights) {
    tryCatch(
      {
        if (layout_type %in% c("stress", "graphopt", "fr")) {
          ggraph::create_layout(graph, layout = layout_type, weights = weights)
        } else {
          ggraph::create_layout(graph, layout = layout_type)
        }
      },
      error = function(e) {
        warning("Layout failed. Falling back to 'fr'.")
        ggraph::create_layout(graph, layout = "fr", weights = weights)
      }
    )
  }
  layout_data <- safely_layout(graph, layout_type, abs(edge_weights))

  # Color handling
  num_colors <- length(unique(communities))
  if (requireNamespace("DspikeIn", quietly = TRUE)) {
    node_colors <- DspikeIn::color_palette$cool_MG[seq_len(num_colors)]
  } else {
    warning("DspikeIn color palette not available. Using rainbow colors.")
    node_colors <- rainbow(num_colors)
  }
  names(node_colors) <- unique(communities)
  igraph::V(graph)$color <- node_colors[as.character(communities)]

  # Label top degree nodes
  top_nodes <- order(-igraph::V(graph)$degree)[seq_len(min(15, igraph::vcount(graph)))]
  label_nodes <- ifelse(igraph::V(graph)$name %in% igraph::V(graph)$name[top_nodes], igraph::V(graph)$name, "")

  # Edge thickness
  edge_thickness <- scales::rescale(abs(edge_weights), to = c(0.5, 3))

  # Generate plot
  plot <- ggraph::ggraph(layout_data) +
    ggraph::geom_edge_link(ggplot2::aes(x = x, y = y, width = edge_thickness), color = "#4F4A4A", alpha = 0.3) +
    ggraph::geom_node_point(ggplot2::aes(x = x, y = y, color = factor(communities), size = igraph::V(graph)$degree), alpha = 0.9, stroke = 1) +
    ggrepel::geom_text_repel(ggplot2::aes(x = x, y = y, label = label_nodes), size = 3, max.overlaps = 15, force = 3) +
    ggplot2::scale_color_manual(values = node_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm")
    ) +
    ggplot2::labs(title = "Modularity-Based Network", color = "Community", size = "Node Degree")

  return(list(
    plot = plot,
    metrics = global_metrics,
    layout_data = layout_data,
    graph = graph
  ))
}


# Usage Example:
# External GraphML file can be loaded using absolute address
# Preloaded internal dataset
# network <- load_graphml("Complete.graphml")
# result <- degree_network(graph_path = network, save_metrics = TRUE, layout_type = "fr")
# print(result$metrics)
# print(result$plot)

# result <- degree_network(graph_path = network , save_metrics = TRUE, layout_type = "circle")
#
