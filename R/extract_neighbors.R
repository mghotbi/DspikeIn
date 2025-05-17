#' @title Extract First and Second Neighbors of a Target Node
#'
#' @description Extracts the first and second neighbors of a given target node in an `igraph` network.
#' Users can provide either an `igraph` object or a GraphML file path (internal or external).
#'
#' @param graph Either:
#'   - An `igraph` object representing a network, OR
#'   - A character string specifying the **file path** to a GraphML file, OR
#'   - `NULL` to use `"Complete.graphml"` from `DspikeIn`.
#' @param target_node Character. The **name** of the target node.
#' @param mode Character. Direction of edges to consider (`"all"`, `"out"`, `"in"`).
#'   - `"all"`: Incoming + outgoing edges (default for undirected graphs).
#'   - `"out"`: Only outgoing edges.
#'   - `"in"`: Only incoming edges.
#'
#' @return A list containing:
#'   \item{first_neighbors}{Character vector of **first-degree** neighbor names.}
#'   \item{second_neighbors}{Character vector of **second-degree** neighbor names.}
#'   \item{summary}{Data frame summarizing the extracted neighbors.}
#'
#' @importFrom igraph neighbors V read_graph
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   Complete <- load_graphml("Complete.graphml")
#'
#'   # Extract neighbors of a specific node from the loaded GraphML file
#'   result1 <- extract_neighbors(graph = Complete, target_node = "OTU69:Basidiobolus_sp")
#'   print(result1$summary)
#'
#'   # for an external GraphML file (enter full address please)
#'   
#'   # Use an already loaded igraph object
#'   Nohubs <- load_graphml("Nohubs.graphml")
#'   result3 <- extract_neighbors(graph = Nohubs, target_node = "OTU1:Lilapila_jurana")
#'   print(result3$summary)
#' }
#'
#' @export
extract_neighbors <- function(graph = NULL, target_node, mode = "all") {
  # ========================
  # Load Graph If Needed
  # ========================
  if (is.null(graph)) {
    message("Using default network: 'Complete.graphml'")
    graph <- load_graphml("Complete.graphml") # Load from package
  } else if (is.character(graph)) {
    # Check if graph is an **external** file or a **package** file
    if (file.exists(graph)) {
      message("Loading external GraphML file: ", graph)
      graph <- igraph::read_graph(graph, format = "graphml")
    } else {
      # Try loading it from `DspikeIn`
      message("Loading GraphML from package: ", graph)
      graph <- load_graphml(graph)
    }
  } else if (!inherits(graph, "igraph")) {
    stop("Error: The 'graph' argument must be an igraph object or a valid GraphML file path.")
  }

  # ========================
  # Validate Target Node
  # ========================
  if (!target_node %in% igraph::V(graph)$name) {
    stop("Error: Target node '", target_node, "' not found in the graph.")
  }

  if (!mode %in% c("all", "out", "in")) {
    stop("Error: Mode must be 'all', 'out', or 'in'.")
  }

  # ========================
  # Find First Neighbors
  # ========================
  first_neighbors <- igraph::neighbors(graph, v = target_node, mode = mode)
  first_neighbors_names <- igraph::V(graph)[first_neighbors]$name

  # ========================
  # Find Second Neighbors
  # ========================
  second_neighbors <- unique(unlist(lapply(first_neighbors_names, function(n) {
    igraph::neighbors(graph, v = n, mode = mode)
  })))

  second_neighbors_names <- igraph::V(graph)[second_neighbors]$name

  # Remove the target node and its first neighbors from second neighbors
  second_neighbors_names <- setdiff(second_neighbors_names, c(first_neighbors_names, target_node))

  # ========================
  # Handle Cases with No Neighbors
  # ========================
  if (length(first_neighbors_names) == 0) {
    warning("Warning: No first neighbors found for '", target_node, "'.")
  }
  if (length(second_neighbors_names) == 0) {
    warning("Warning: No second neighbors found for '", target_node, "'.")
  }

  # ========================
  # Create Summary Data Frame
  # ========================
  summary_table <- data.frame(
    Type = c(
      rep("First Neighbor", length(first_neighbors_names)),
      rep("Second Neighbor", length(second_neighbors_names))
    ),
    Node = c(first_neighbors_names, second_neighbors_names),
    stringsAsFactors = FALSE
  )

  # ========================
  # Return Output
  # ========================
  return(list(
    first_neighbors = first_neighbors_names,
    second_neighbors = second_neighbors_names,
    summary = summary_table
  ))
}


# Usage Example:

# Complete <- load_graphml("Complete.graphml")
# result2 <- extract_neighbors(graph = Complete,
# target_node = "OTU69:Basidiobolus_sp", mode = "all")
# print(result3$summary)

# result <- extract_neighbors(graph = network,
# target_node = "OTU1:Lilapila_jurana", mode = "all")

# result <- extract_neighbors(graph = "~/herp.spiecsym.network.graphml",
# target_node = "OTU1:Lilapila_jurana", mode = "all")
