#' @title Compute and Visualize Node-Level Network Metrics
#' @description Computes various network metrics and generates:
#' - Two **multi-panel visualizations** (each with 4 subplots).
#' - A **faceted plot** of standardized (Z-score) metrics across communities.
#' - A **formatted flextable** summarizing node-level metrics.
#'
#' @details The function computes the following node-level metrics:
#'
#' | **Metric**                     | **Description** |
#' |--------------------------------|-----------------------------------------------------------|
#' | `Node`                         | Node name (character format) |
#' | `Degree`                       | Number of edges connected to the node |
#' | `Strength`                     | Sum of edge weights connected to the node |
#' | `Closeness`                    | Closeness centrality (normalized, based on shortest paths) |
#' | `Betweenness`                  | Betweenness centrality (normalized, measures control over network flow) |
#' | `EigenvectorCentrality`        | Eigenvector centrality (importance based on connections to influential nodes) |
#' | `PageRank`                     | PageRank score (importance based on incoming links) |
#' | `Transitivity`                 | Local clustering coefficient (tendency of a node to form triangles) |
#' | `Coreness`                     | Node's coreness (from k-core decomposition) |
#' | `Constraint`                   | Burt's constraint (measures structural holes in a node's ego network) |
#' | `EffectiveSize`                | Inverse of constraint (larger values = more non-redundant connections) |
#' | `Redundancy`                   | Sum of constraint values of a node's alters |
#' | `Community`                    | Community assignment from Louvain clustering |
#' | `Efficiency`                   | Global efficiency (average inverse shortest path length) |
#' | `Local_Efficiency`             | Local efficiency (subgraph efficiency for a node's neighbors) |
#' | `Within_Module_Connectivity`   | Proportion of neighbors in the same community |
#' | `Among_Module_Connectivity`    | Proportion of neighbors in different communities |
#'
#' @param graph An `igraph` object representing the network.
#' @param save_path Character. A file path (without extension) to save figures and table. Default is NULL (no saving).
#'
#' @return A list containing:
#' - `metrics`: A data frame with node-level metrics.
#' - `plot1`: First multi-panel plot (2x2 layout with 4 subplots)
#' - `plot2`: Second multi-panel plot (2x2 layout with 4 subplots)
#' - `facet_plot`: A faceted plot showing **Z-score standardized** metrics across communities.
#'
#' @importFrom igraph V E degree strength closeness betweenness eigen_centrality page_rank transitivity coreness
#' @importFrom igraph cluster_louvain constraint neighbors induced_subgraph distances
#' @importFrom ggplot2 ggplot aes geom_histogram geom_point geom_boxplot labs theme_minimal scale_color_viridis_c
#' @importFrom ggplot2 scale_fill_viridis_d ggsave scale_color_manual facet_grid element_blank theme
#' @importFrom patchwork wrap_plots
#' @importFrom reshape2 melt
#' @importFrom dplyr mutate_at
#' @examples
#' library(igraph)
#' set.seed(42)
#' # For external graphml please use full address
#' # Load internal graphml
#' Complete <- load_graphml("Complete.graphml")
#'
#' # Compute node-level metrics
#' result <- node_level_metrics(Complete)
#'
#' # View computed metrics
#' print(result$metrics)
#'
#' # Show the first 4x4 plot
#' print(result$plots$plot1)
#'
#' # Show the second 4x4 plot
#' print(result$plots$plot2)
#'
#' # Show facet plot
#' print(result$facet_plot)
#'
#' # Print metrics and flextable
#' print(result$metrics)
#' print(result$flextable)
#' @export
node_level_metrics <- function(graph, save_path = NULL) {
  # =====================
  # Ensure Edge Weights
  # =====================
  if (!"weight" %in% names(igraph::edge.attributes(graph))) {
    igraph::E(graph)$weight <- rep(1, igraph::ecount(graph))
  } else {
    igraph::E(graph)$weight <- abs(igraph::E(graph)$weight)
  }

  edge_weights <- igraph::E(graph)$weight
  clusters <- igraph::cluster_louvain(graph, weights = edge_weights)

  # =====================
  # Compute Node-Level Metrics
  # =====================
  metrics <- data.frame(
    Node = as.character(igraph::V(graph)$name),
    Degree = as.integer(igraph::degree(graph)),
    Strength = igraph::strength(graph, weights = edge_weights),
    Closeness = igraph::closeness(graph, weights = edge_weights, normalized = TRUE),
    Betweenness = withCallingHandlers(
      igraph::betweenness(graph, weights = edge_weights, normalized = TRUE),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    EigenvectorCentrality = igraph::eigen_centrality(graph, weights = edge_weights)$vector,
    PageRank = igraph::page_rank(graph, weights = edge_weights)$vector,
    Harmonic = igraph::harmonic_centrality(graph, weights = edge_weights),
    Transitivity = igraph::transitivity(graph, type = "local", isolates = "zero"),
    Coreness = as.integer(igraph::coreness(graph, mode = "all")),
    Constraint = igraph::constraint(graph),
    EffectiveSize = 1 / igraph::constraint(graph),
    Redundancy = vapply(igraph::V(graph), function(v) {
      sum(igraph::constraint(graph)[v])
    }, numeric(1)),
    Community = as.factor(clusters$membership),
    Efficiency = NA,
    Local_Efficiency = NA,
    Within_Module_Connectivity = NA,
    Among_Module_Connectivity = NA
  )

  # =====================
  # Compute Efficiency & Local Efficiency
  # =====================
  metrics$Efficiency <- vapply(igraph::V(graph), function(v) {
    dist_mat <- igraph::distances(graph, v = v)
    valid_distances <- dist_mat[dist_mat > 0 & is.finite(dist_mat)]
    if (length(valid_distances) == 0) {
      return(0)
    }
    mean(1 / valid_distances, na.rm = TRUE)
  }, numeric(1))

  metrics$Local_Efficiency <- vapply(igraph::V(graph), function(v) {
    neighbors <- igraph::neighbors(graph, v)
    if (length(neighbors) < 2) {
      return(0)
    }
    subgraph <- igraph::induced_subgraph(graph, neighbors)
    dist_mat <- igraph::distances(subgraph)
    valid_distances <- dist_mat[dist_mat > 0 & is.finite(dist_mat)]
    if (length(valid_distances) == 0) {
      return(0)
    }
    mean(1 / valid_distances, na.rm = TRUE)
  }, numeric(1))

  # =====================
  # Within/Among Module Connectivity
  # =====================
  metrics$Within_Module_Connectivity <- vapply(igraph::V(graph), function(v) {
    node_community <- clusters$membership[as.numeric(v)]
    neighbors <- as.numeric(igraph::neighbors(graph, v))
    if (length(neighbors) == 0) {
      return(0)
    }
    same_community_neighbors <- sum(clusters$membership[neighbors] == node_community)
    same_community_neighbors / length(neighbors)
  }, numeric(1))

  metrics$Among_Module_Connectivity <- vapply(igraph::V(graph), function(v) {
    node_community <- clusters$membership[as.numeric(v)]
    neighbors <- as.numeric(igraph::neighbors(graph, v))
    if (length(neighbors) == 0) {
      return(0)
    }
    diff_community_neighbors <- sum(clusters$membership[neighbors] != node_community)
    diff_community_neighbors / length(neighbors)
  }, numeric(1))

  # =====================
  # Generate Plots
  # =====================
  color_palette <- DspikeIn::color_palette$cool_MG

  plot1 <- patchwork::wrap_plots(
    ggplot2::ggplot(metrics, ggplot2::aes(x = Degree)) +
      ggplot2::geom_histogram(fill = "#FB5607", bins = 20, alpha = 0.7) +
      ggplot2::theme_minimal(),
    ggplot2::ggplot(metrics, ggplot2::aes(x = PageRank)) +
      ggplot2::geom_histogram(fill = "#553C9A", bins = 20, alpha = 0.7) +
      ggplot2::theme_minimal(),
    ggplot2::ggplot(metrics, ggplot2::aes(y = Closeness, x = Betweenness, color = Community)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top"),
    ggplot2::ggplot(metrics, ggplot2::aes(y = Strength, x = EigenvectorCentrality, color = Degree)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme_minimal(),
    ncol = 2
  )

  plot2 <- patchwork::wrap_plots(
    ggplot2::ggplot(metrics, ggplot2::aes(y = Among_Module_Connectivity, x = Betweenness, color = Community)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top"),
    ggplot2::ggplot(metrics, ggplot2::aes(y = Redundancy, x = Degree, color = Community)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top"),
    ggplot2::ggplot(metrics, ggplot2::aes(y = Constraint, x = Within_Module_Connectivity, color = Community)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top"),
    ggplot2::ggplot(metrics, ggplot2::aes(y = Degree, x = Coreness, color = Community)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::scale_color_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top"),
    ncol = 2
  )

  # =====================
  # Generate Facet Plot
  # =====================
  metrics_zscore <- dplyr::mutate_if(metrics, is.numeric, as.double)
  metrics_zscore <- dplyr::mutate_at(metrics_zscore, dplyr::vars(-Node, -Community), ~ scale(.))
  metrics_long <- reshape2::melt(metrics_zscore,
    id.vars = c("Node", "Community"),
    variable.name = "Metric", value.name = "Z_Score"
  )

  facet_plot <- ggplot2::ggplot(metrics_long, ggplot2::aes(x = Community, y = Z_Score, fill = Community)) +
    ggplot2::geom_boxplot(alpha = 0.8) +
    ggplot2::facet_grid(~Metric, scales = "free_y") +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  # =====================
  # Generate Flextable
  # =====================
  unique_communities <- levels(metrics$Community)
  community_colors <- setNames(color_palette[seq_len(length(unique_communities))], unique_communities)
  bg_colors <- community_colors[as.character(metrics$Community)]

  table_flex <- flextable::qflextable(metrics) |>
    flextable::theme_booktabs() |>
    flextable::set_table_properties(layout = "autofit") |>
    flextable::bold(part = "header") |>
    flextable::color(j = "Community", color = "white", part = "body") |>
    flextable::bg(j = "Community", bg = bg_colors, part = "body") |>
    flextable::autofit()

  return(list(
    metrics = metrics,
    flextable = table_flex,
    plot1 = plot1,
    plot2 = plot2,
    facet_plot = facet_plot
  ))
}


# Usage Example:
# set.seed(42)
# CustomNet <- load_graphml("~/herp.spiecsym.network.graphml")
# result <- node_level_metrics(CustomNet)
# result$facet_plot

# load internal network
# Complete <- load_graphml("Complete.graphml")
# result <- node_level_metrics(Complete)
# result$metrics
# View computed metrics
# print(result$metrics)

# Show the first 4x4 plot
# print(result$plots$plot1)

# Show the second 4x4 plot
# print(result$plots$plot2)

# result$table
# ggplot2::ggsave("plot1.png", result$plot1, width = 10, height = 5)
# ggplot2::ggsave("plot2.png", result$plot2, width = 10, height = 5)
# ggplot2::ggsave("facet_plot.png", result$facet_plot, width = 12, height = 6)
# 
