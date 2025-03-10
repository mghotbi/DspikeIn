#' @title Analyze and Visualize a Microbial Network
#'
#' @description This function loads a microbial network from a GraphML file or an `igraph` object,
#' computes node degree and modularity, assigns node sizes based on degree, and visualizes the network.
#' Edge thickness is determined by weight. It also computes and saves global network metrics.
#'
#' @param graph_path Character or igraph object. Path to the GraphML file containing the network,
#' or an already loaded `igraph` object.
#' @param save_metrics Logical. If TRUE, saves the global metrics as `"Global_Network_Metrics.csv"`.
#' @param layout_type Character. The layout algorithm for visualization. Options: `"stress"` (default), `"graphopt"`,
#' `"fr"` (Fruchterman-Reingold), `"mds"`, `"kk"` (Kamada-Kawai).
#'
#' @return A list containing:
#'   \item{plot}{A ggplot object displaying the network.}
#'   \item{metrics}{A dataframe with global network metrics.}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   Complete <- load_graphml("Complete.graphml")
#'
#'   # Compute degree metrics and visualize the network
#'   result <- degree_network(graph_path = Complete, save_metrics = TRUE)
#'   print(result$metrics)
#'   print(result$plot)
#'
#'   # Use Kamada-Kawai layout for better separation
#'   # Options: `"stress"` (default), `"graphopt"`, `"fr"`
#'   # (Fruchterman-Reingold), `"mds"`, `"kk"` (Kamada-Kawai)
#'   result_kk <- degree_network(
#'     graph_path = load_graphml("Complete.graphml"),
#'     save_metrics = TRUE,
#'     layout_type = "kk"
#'   )
#'   print(result_kk$plot)
#' }
#' }
#' @importFrom igraph read_graph cluster_louvain membership degree edge_density is_connected components induced_subgraph V E
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales rescale
#' @importFrom ggplot2 theme_void labs theme scale_color_manual element_rect guides guide_legend
#' @importFrom utils write.csv
#' @export
degree_network <- function(graph_path, save_metrics = TRUE, layout_type = "stress") {

  # ================
  #  Load Network Graph
  # ================
  if (inherits(graph_path, "igraph")) {
    message("\U0001F504 Using preloaded igraph object...")
    graph <- graph_path
  } else {
    if (!file.exists(graph_path)) stop("\U0000274C Error: Graph file not found.")
    message("\U0001F504 Loading graph from file...")
    graph <- igraph::read_graph(graph_path, format = "graphml")
  }

  # ================
  #  Process Edge Weights
  # ================
  message("\U0001F503 Processing edge weights...")
  if (!is.null(igraph::E(graph)$weight)) {
    igraph::E(graph)$original_weight <- igraph::E(graph)$weight  # Store original values
    igraph::E(graph)$weight <- ifelse(igraph::E(graph)$weight <= 0, 1e-3, igraph::E(graph)$weight)  # Ensure positive weights
  } else {
    igraph::E(graph)$weight <- rep(1, igraph::ecount(graph))  # Assign default weight if missing
  }

  # ================
  #  Compute Node Degree
  # ================
  message("\U0001F504 Computing node degree...")
  igraph::V(graph)$degree <- igraph::degree(graph)  # Store node degree

  # ================
  #  Compute Network Modularity & Metrics
  # ================
  message("\U0001F310 Computing global network metrics...")
  edge_weights <- igraph::E(graph)$original_weight
  is_conn <- igraph::is_connected(graph)

  largest_comp <- if (!is_conn) {
    igraph::induced_subgraph(graph, which(igraph::components(graph)$membership == which.max(igraph::components(graph)$csize)))
  } else {
    graph
  }

  global_metrics <- data.frame(
    Nodes = igraph::vcount(graph),
    Edges = igraph::ecount(graph),
    Modularity = tryCatch(igraph::modularity(igraph::cluster_louvain(graph)), error = function(e) NA),
    Density = tryCatch(igraph::edge_density(graph), error = function(e) NA),
    Transitivity = tryCatch(igraph::transitivity(graph), error = function(e) NA),
    Diameter = tryCatch(if (is_conn) igraph::diameter(largest_comp, weights = edge_weights) else NA, error = function(e) NA),
    Avg_Path_Length = tryCatch(if (is_conn) igraph::mean_distance(largest_comp, weights = edge_weights) else NA, error = function(e) NA),
    Avg_Degree = tryCatch(mean(igraph::degree(graph)), error = function(e) NA)
  )

  # Save global metrics as CSV if required
  if (save_metrics) {
    write.csv(global_metrics, "Global_Network_Metrics.csv", row.names = FALSE)
    message("\U0001F4C2 Global network metrics saved as 'Global_Network_Metrics.csv'.")
  }

  # ================
  #  Validate Layout & Generate Layout
  # ================
  valid_layouts <- c("stress", "graphopt", "fr", "mds", "kk")
  if (!layout_type %in% valid_layouts) {
    stop("Invalid layout. Choose from: 'stress', 'graphopt', 'fr', 'mds', 'kk'")
  }

  message("\U0001F5A7 Generating network layout using '", layout_type, "'...")
  layout_data <- ggraph::create_layout(graph, layout = layout_type, weights = abs(edge_weights))

  # ================
  #  Assign Colors and Labels
  # ================
  edge_thickness <- scales::rescale(abs(edge_weights), to = c(0.5, 3))  # Scaled edge thickness based on weight

  # Node Colors Based on Community
  communities <- tryCatch(igraph::membership(igraph::cluster_louvain(graph)), error = function(e) rep(1, igraph::vcount(graph)))

  # Use the Internal Color Palette from DspikeIn
  num_colors <- length(unique(communities))
  node_colors <- tryCatch(DspikeIn::color_palette$MG[1:num_colors], error = function(e) rainbow(num_colors))

  # Ensure colors match communities
  names(node_colors) <- unique(communities)
  igraph::V(graph)$color <- node_colors[as.character(communities)]

  # Label Top 15 Nodes by Degree
  top_nodes <- order(-igraph::V(graph)$degree)[1:min(15, vcount(graph))]
  label_nodes <- ifelse(seq_along(igraph::V(graph)$name) %in% top_nodes, igraph::V(graph)$name, "")

  # ================
  #  Generate Plot
  # ================
  plot <- ggraph::ggraph(layout_data) +
    ggraph::geom_edge_link(ggplot2::aes(x = x, y = y, width = edge_thickness), color =  "#4F4A4A" , alpha = 0.3) +
    ggraph::geom_node_point(ggplot2::aes(x = x, y = y, color = factor(communities), size = igraph::V(graph)$degree), alpha = 0.9, stroke = 1) +
    ggrepel::geom_text_repel(ggplot2::aes(x = x, y = y, label = label_nodes), size = 3, max.overlaps = 15, force = 3) +
    ggplot2::scale_color_manual(values = node_colors) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom", legend.key.size = ggplot2::unit(0.4, "cm")) +
    ggplot2::labs(title = "Modularity-Based Network", color = "Community", size = "Node Degree")

  return(list(plot = plot, metrics = global_metrics))
}

# Usage Example:
# External GraphML file can be loaded using absolute address
# Preloaded internal dataset
# network <- load_graphml("Complete.graphml")
# result <- degree_network(graph_path = network, save_metrics = TRUE, layout_type = "fr")
# print(result$metrics)
# print(result$plot)

# result <- degree_network(graph_path = network , save_metrics = TRUE, layout_type = "circle")
