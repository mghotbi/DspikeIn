#' @title Simulate Network Robustness under Node Removal
#'
#' @description Evaluates the robustness of a network by iteratively removing nodes
#' and measuring the size of the largest connected component at each step.
#'
#' @details Users can specify a node removal strategy from the following options:
#'   - **"random"**: Removes nodes randomly.
#'   - **"degree"**: Removes the node with the highest degree first.
#'   - **"betweenness"**: Removes the node with the highest betweenness centrality first.
#'
#' **Important:** For the `"betweenness"` strategy, all edge weights **must be positive**.
#'
#' @param graph An `igraph` object representing the network.
#' @param steps Integer. Number of nodes to remove. Default is `10`.
#' @param removal_strategy Character. Node removal strategy, one of:
#'   - `"random"` (default): Removes nodes randomly.
#'   - `"degree"`: Removes the highest-degree nodes first.
#'   - `"betweenness"`: Removes the highest-betweenness nodes first.
#' @param plot_results Logical. If `TRUE`, generates and **returns** a plot of the robustness curve. Default is `TRUE`.
#'
#' @return A **list** containing:
#'   - `results`: A `data.frame` with the step number and the relative size of the largest connected component.
#'   - `plot`: A `ggplot2` object (returned if `plot_results = TRUE`).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   Complete <- load_graphml("Complete.graphml")
#'
#'   # Simulate robustness by removing 200 highest-degree nodes
#'   robustness_degree <- simulate_network_robustness(
#'     graph = Complete,
#'     steps = 200,
#'     removal_strategy = "degree"
#'   )
#'
#'   # Simulate robustness with random node removal
#'   robustness_random <- simulate_network_robustness(
#'     graph = Complete,
#'     steps = 200,
#'     removal_strategy = "random"
#'   )
#'
#'   # Simulate robustness with betweenness-based node removal
#'   robustness_betweenness <- simulate_network_robustness(
#'     graph = Complete,
#'     steps = 200,
#'     removal_strategy = "betweenness"
#'   )
#'
#'   # Print robustness plots
#'   print(robustness_degree$plot)
#'   print(robustness_random$plot)
#'   print(robustness_betweenness$plot)
#' }
#' }
#' @importFrom igraph vcount ecount components delete_vertices degree betweenness V E is_weighted
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal
#'
#' @export
simulate_network_robustness <- function(graph, steps = 10, removal_strategy = "random", plot_results = TRUE) {
  message("\U0001F578 Running network robustness simulation...")

  # =====================
  # Validate Input Graph
  # =====================
  if (!inherits(graph, "igraph")) {
    stop("\U0000274C Error: Input is not a valid igraph object.")
  }
  if (igraph::vcount(graph) == 0 || igraph::ecount(graph) == 0) {
    stop("\U0000274C Error: The graph is empty or has no edges.")
  }


  # =====================
  # Ensure Positive Edge Weights for Betweenness Calculation
  # =====================
  if (removal_strategy == "betweenness" && igraph::is_weighted(graph)) {
    if (any(igraph::E(graph)$weight <= 0)) {
      warning("\U00026A0 Adjusting edge weights to be positive for betweenness calculation.")
      igraph::E(graph)$weight <- igraph::E(graph)$weight + abs(min(igraph::E(graph)$weight)) + 1e-5
    }
  }

  # =====================
  # Initialize Variables
  # =====================
  g <- graph  # Copy the graph to avoid modifying the original
  results <- data.frame(Step = integer(steps), LargestComponentRatio = numeric(steps))

  # =====================
  # Node Removal Simulation
  # =====================
  for (i in seq_len(steps)) {
    if (igraph::vcount(g) == 0) break  # Stop if the graph is empty

    # Compute the largest component size as a fraction of original nodes
    largest_component <- max(igraph::components(g)$csize)
    results[i, ] <- c(i, largest_component / igraph::vcount(graph))

    # Determine node removal strategy
    if (removal_strategy == "random") {
      remove_node <- sample(igraph::V(g), 1)
    } else if (removal_strategy == "degree") {
      remove_node <- igraph::V(g)[which.max(igraph::degree(g))]
    } else if (removal_strategy == "betweenness") {
      if (igraph::is_weighted(g)) {
        betweenness_values <- igraph::betweenness(g, weights = igraph::E(g)$weight)
      } else {
        betweenness_values <- igraph::betweenness(g)
      }
      remove_node <- igraph::V(g)[which.max(betweenness_values)]
    } else {
      stop("\U0000274C Error: Invalid removal strategy. Choose 'random', 'degree', or 'betweenness'.")
    }

    # Remove the selected node
    g <- igraph::delete_vertices(g, remove_node)
  }

  # =====================
  # Generate Robustness Plot (Optional)
  # =====================
  robustness_plot <- NULL
  if (plot_results) {
    message("\U0001F504 Generating robustness plot...")
    robustness_plot <- ggplot2::ggplot(results, ggplot2::aes(x = Step, y = LargestComponentRatio)) +
      ggplot2::geom_line(color = "blue4", linewidth = 1) +
      ggplot2::labs(
        title = "Network Robustness Simulation",
        x = "Step (Nodes Removed)",
        y = "Largest Component Ratio"
      ) +
      ggplot2::theme_minimal()
  }

  return(list(results = results, plot = robustness_plot))
}


# =====================
# Example Usage
# =====================
# Load a network from a GraphML file
# Complete <- load_graphml("Complete.graphml")
# Simulate robustness by removing 200 highest-degree nodes
# robustness_results <- simulate_network_robustness(graph = NoBasid,
# steps = 200, removal_strategy = "degree")
 # Simulate robustness by removing 200 highest-degree nodes
# robustness_results <- simulate_network_robustness(graph = NoHubs,
# steps = 200, removal_strategy = "degree")
# robustness_results <- simulate_network_robustness(graph = Complete,
# steps = 200, removal_strategy = "betweenness")

# robustness_results <- simulate_network_robustness(graph = Complete,
# steps = 200, removal_strategy = "betweenness")

# print(head(robustness_results$results))

# Display plot
# print(robustness_results$plot)


