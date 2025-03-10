#' @title Plot a Phylogenetic Tree with Glommed OTUs
#' @description This function gloms the OTUs in a `phyloseq` object at the specified resolution and plots the phylogenetic tree.
#' If a `TreeSummarizedExperiment` (TSE) object is provided, it is first converted to `phyloseq`.
#'
#' @param physeq A `phyloseq` or `TreeSummarizedExperiment` object containing the phylogenetic data.
#' @param resolution A numeric value specifying the resolution for glomming the OTUs.
#' @param output_prefix A character string specifying the prefix for the output files. Default is `"glommed_tree"`.
#' @param width A numeric value specifying the width of the output plot. Default is `20`.
#' @param height A numeric value specifying the height of the output plot. Default is `18`.
#' @param print_plot A logical value specifying whether to print the plot. Default is `TRUE`.
#' @return `NULL`. The function saves the plot to the specified output file.
#'
#' @importFrom phyloseq phy_tree
#' @importFrom speedyseq tree_glom
#' @importFrom ggtree ggtree geom_tiplab theme_tree2 geom_text2 geom_nodepoint
#' @importFrom ggplot2 ggsave aes ggtitle
#' @importFrom ape Nnode Ntip
#' @importFrom dplyr filter mutate
#' @importFrom tidytree as_tibble
#' @examples
#' \dontrun{
#' plot_glommed_tree(Tetragenococcus_tse, resolution = 0.2, output_prefix = "top",
#'  width = 18, height = 18, print_plot = TRUE)
#'
#' plot_glommed_tree(Tetragenococcus, resolution = 0.3, output_prefix = "tse_tree")
#' }
#' @export
plot_glommed_tree <- function(physeq, resolution, output_prefix = "glommed_tree", width = 20, height = 18, print_plot = TRUE) {

  # Convert TSE to Phyloseq if needed
  if (inherits(physeq, "TreeSummarizedExperiment")) {
    message(" Detected TSE object. Converting to phyloseq for tree glomming...")
    physeq <- convert_tse_to_phyloseq(physeq)
  }

  # Perform tree glomming
  physeq_glommed <- speedyseq::tree_glom(physeq, resolution = resolution)

  # Check if tree glomming worked correctly
  if (is.null(physeq_glommed)) {
    stop(" Tree glomming failed. Please check the input phyloseq object and resolution.")
  }

  # Extract tree from phyloseq object
  phy_tree <- phyloseq::phy_tree(physeq_glommed)

  # Convert tree structure for ggtree
  tree_data <- ggtree::fortify(phy_tree)

  # 🛠 Ensure bootstrap labels are only assigned to internal nodes
  num_tips <- ape::Ntip(phy_tree)  # Number of tips
  num_nodes <- ape::Nnode(phy_tree)  # Number of internal nodes
  node_labels <- phy_tree$node.label  # Bootstrap values

  # Filter valid bootstrap labels (remove NA)
  if (!is.null(node_labels)) {
    node_labels <- ifelse(node_labels == "", NA, node_labels)  # Convert empty strings to NA
    valid_labels <- !is.na(node_labels)

    if (sum(valid_labels) > 0) {
      message(" Adding bootstrap values to internal nodes.")

      # Get only internal nodes (node numbers start after num_tips)
      bootstrap_df <- tree_data %>%
        dplyr::filter(node > num_tips) %>%
        dplyr::mutate(label = node_labels[valid_labels])  # Assign correct bootstrap values

      # Create ggtree object
      tree_plot <- ggtree::ggtree(phy_tree) +
        ggtree::geom_tiplab(ggplot2::aes(label = label), size = 3) +  # Label tips
        ggtree::theme_tree2() +
        ggplot2::ggtitle(paste("Glommed Tree at Resolution", resolution)) +
        ggtree::geom_nodepoint(size = 3, color = "blue")  # Highlight internal nodes

      # Add bootstrap labels only to internal nodes with correct `x` and `y` mapping
      tree_plot <- tree_plot +
        ggtree::geom_text2(data = bootstrap_df, ggplot2::aes(x = x, y = y, label = label),
                           size = 3, vjust = -0.5, color = "red4")

    } else {
      warning(" No valid bootstrap values found. Skipping bootstrap labels.")

      # Create ggtree object without bootstrap labels
      tree_plot <- ggtree::ggtree(phy_tree) +
        ggtree::geom_tiplab(ggplot2::aes(label = label), size = 3) +  # Label tips
        ggtree::theme_tree2() +
        ggplot2::ggtitle(paste("Glommed Tree at Resolution", resolution)) +
        ggtree::geom_nodepoint(size = 3, color = "blue")  # Highlight internal nodes
    }
  }

  # Print the plot if required
  if (print_plot) {
    print(tree_plot)
  }

  # Save the plot
  ggplot2::ggsave(paste0(output_prefix, "_tree.pdf"), plot = tree_plot, width = width, height = height)
  cat(" Glommed tree plot saved to:", paste0(output_prefix, "_tree.pdf"), "\n")

  return(NULL)
}


# Example usage:
#plot_glommed_tree(Tetragenococcus_tse, resolution = 0.2, output_prefix = "top", width = 18, height = 18, print_plot = TRUE)
#plot_glommed_tree(Tetragenococcus, resolution = 0.3, output_prefix = "tser_tree")
