#' @title Plot a Phylogenetic Tree with Multiple Sequence Alignment
#'
#' @description Generates a phylogenetic tree plot with an additional multiple sequence alignment (MSA) visualization.
#' Supports both `phyloseq` and `TreeSummarizedExperiment` (TSE) objects.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing the phylogenetic data.
#' @param output_prefix A character string specifying the prefix for the output files. Default is `"physeq_tree_alignment"`.
#' @param width A numeric value specifying the width of the output plot. Default is `15`.
#' @param height A numeric value specifying the height of the output plot. Default is `15`.
#' @return `NULL`. The function saves the plots as PNG and PDF files.
#'
#' @details
#' - Extracts the **phylogenetic tree** and **reference sequences** from the input object.
#' - Performs **multiple sequence alignment (MSA)** on the reference sequences.
#' - Plots the tree and **aligns it with the MSA** using `ggtree::msaplot()`.
#' - Saves the output as **PNG and PDF files**.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Generate and save a tree plot with alignment
#'   plot_tree_with_alignment(
#'     physeq_16SOTU,
#'     output_prefix = "tree_alignment",
#'     width = 15,
#'     height = 15
#'   )
#'
#'   # Example usage with a TreeSummarizedExperiment (TSE) object:
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   plot_tree_with_alignment(
#'     tse_16SOTU,
#'     output_prefix = "tse_tree_alignment",
#'     width = 15,
#'     height = 15
#'   )
#' }
#' }
#'
#' @importFrom phyloseq phy_tree refseq
#' @importFrom DECIPHER AlignSeqs
#' @importFrom Biostrings writeXStringSet
#' @importFrom ggtree ggtree theme_tree2 msaplot geom_tiplab geom_tippoint
#' @importFrom ggplot2 scale_color_continuous theme ggsave
#' @importFrom grid unit
#' @importFrom S4Vectors metadata
#' @export
plot_tree_with_alignment <- function(obj, output_prefix = "physeq_tree_alignment", width = 15, height = 15) {
  suppressMessages({
    if (!inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
      stop("\U0000274C Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
    }

    # Extract phylogenetic tree
    tree <- get_phy_tree(obj)
    if (is.null(tree) || is.null(tree$edge.length) || length(tree$edge.length) == 0) {
      stop("Tree has no valid branch lengths.")
    }

    # Extract reference sequences
    ref_sequences <- get_reference_seq(obj)
    if (is.null(ref_sequences)) {
      stop("Reference sequences not found in the object.")
    }

    # Perform multiple sequence alignment (MSA)
    alignment <- DECIPHER::AlignSeqs(ref_sequences, anchor = NA)

    # Write aligned sequences to a temporary FASTA file
    temp_fasta <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(alignment, temp_fasta)

    # Generate ggtree plot
    tree_plot <- ggtree::ggtree(tree, layout = "fan", branch.length = "branch.length", ladderize = TRUE,
                                right = TRUE, root.position = 0) + ggtree::theme_tree2()

    # Add multiple sequence alignment plot
    ms <- ggtree::msaplot(p = tree_plot, fasta = temp_fasta)

    # Customize the plot appearance
    final_plot <- ms + ggplot2::scale_color_continuous(low = 'gray88', high = 'red4') +
      ggtree::geom_tiplab(ggplot2::aes(color = branch.length), align = TRUE, size = 3) +
      ggplot2::theme(legend.position = "right", legend.key.size = grid::unit(3, "lines")) +
      ggtree::geom_tippoint(color = "#BB650B", shape = "*", size = 3) + ggtree::theme_tree2()

    print(final_plot)

    # Save the plots
    png_filename <- paste0(output_prefix, ".png")
    pdf_filename <- paste0(output_prefix, ".pdf")

    ggplot2::ggsave(png_filename, plot = final_plot, width = width, height = height, units = "in")
    ggplot2::ggsave(pdf_filename, plot = final_plot, width = width, height = height, units = "in")

    cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
  })
}

# Example usage with a phyloseq object
# plot_tree_with_alignment(Tet, output_prefix = "tree_alignment", width = 15, height = 15)

# Example usage with a TreeSummarizedExperiment object
#plot_tree_with_alignment(tse, output_prefix = "tse_tree_alignment", width = 15, height = 15)
