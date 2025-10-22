#' @title Load GraphML Without 'id' Conflicts
#' @description Loads a GraphML network from `DspikeIn` or user path, ensuring all attributes remain intact and correctly assigning original node names.
#'
#' @param filename Name of the GraphML file (default: "herp.spiecsym.network.graphml") or full file path.
#' @return An `igraph` object representing the network.
#'
#' @importFrom igraph read_graph set_vertex_attr delete_vertex_attr vertex_attr_names V
#' @importFrom xml2 read_xml write_xml xml_find_all xml_attr xml_set_attr xml_remove
#' @export
load_graphml <- function(filename = "herp.spiecsym.network.graphml") {
  # If the file exists externally, use it directly
  if (file.exists(filename)) {
    file_path <- filename
    message("Loading external GraphML file: ", filename)
  } else {
    file_path <- system.file("extdata", filename, package = "DspikeIn")

    if (file_path == "") {
      stop("Error: The specified GraphML file does not exist in the package or as a user file.")
    }

    message("Loading package GraphML file: ", filename)
  }
  # Suppress known attribute warning on import
  graph <- withCallingHandlers(
    igraph::read_graph(file_path, format = "graphml"),
    warning = function(w) invokeRestart("muffleWarning")
  )
  # Step 1: Identify the correct attribute for node names
  name_attr <- NULL
  if ("label" %in% igraph::vertex_attr_names(graph)) {
    name_attr <- "label"
  } else if ("name" %in% igraph::vertex_attr_names(graph)) {
    name_attr <- "name"
  } else if ("id" %in% igraph::vertex_attr_names(graph)) {
    name_attr <- "id"
  }

  # Step 2: Assign correct node names
  if (!is.null(name_attr)) {
    graph <- igraph::set_vertex_attr(graph, "name", value = igraph::vertex_attr(graph, name_attr))
    message("Assigned correct node names from attribute: ", name_attr)
  } else {
    message("Warning: No valid node name attribute found. Keeping default names.")
  }

  # Ensure the graph loaded correctly
  message(
    "Successfully loaded: ", filename,
    " with ", igraph::vcount(graph), " nodes and ", igraph::ecount(graph), " edges."
  )

  return(graph)
}


# Usage Example
# library(DspikeIn)

# Load the internal GraphML files
# Complete <- load_graphml("Complete.graphml")
# Customnet <- load_graphml("~/Customnet.graphml")

# NoHubs <- load_graphml("NoHubs.graphml")
# NoBasid <- load_graphml("NoBasid.graphml")

# Check the graphs
# print(Complete)
# print(NoHubs)
# print(NoBasid)
#
