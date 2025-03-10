#' @title Plot a Neighbor-Joining Phylogenetic Tree with Bootstrap Values
#'
#' @description Constructs a **Neighbor-Joining (NJ) tree** with **bootstrap values**,
#' extracted from a `phyloseq` or `TreeSummarizedExperiment` object, or directly from a FASTA file.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object containing reference sequences,
#' or a **FASTA file** path (`character`) containing DNA sequences.
#' @param output_file A character string specifying the output file for the tree plot.
#' Default is `"neighbor_joining_tree.png"`.
#' @return `NULL`. The function saves the tree plot with bootstrap values to the specified output file.
#'
#' @details
#' - Extracts **reference sequences** from the object (or reads from a FASTA file).
#' - Performs **multiple sequence alignment (MSA)** for proper tree reconstruction.
#' - Constructs a **Neighbor-Joining (NJ) tree** with **bootstrap support values**.
#' - Uses **polished aesthetics** for a **publication-ready tree plot**.
#'
#' @examples
#' \dontrun{
#' data("physeq_16SOTU", package = "DspikeIn")
#' plot_tree_nj(physeq_16SOTU, output_file = "nj_tree_physeq.png")
#'
#' tse_obj <-convert_phyloseq_to_tse(physeq_16SOTU)
#' plot_tree_nj(tse_obj, output_file = "nj_tree_tse.png")
#'
#' # Example usage with your *FASTA* file:
#' sequences <- load_fasta("combined.fasta")
#' plot_tree_nj("sequences", output_file = "nj_tree_fasta.png")
#' }
#'
#' @importFrom ape nj plot.phylo boot.phylo nodelabels ladderize
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom msa msa
#' @importFrom phangorn phyDat dist.ml
#' @importFrom S4Vectors metadata
#' @importFrom grDevices png dev.off
#' @export
plot_tree_nj <- function(obj, output_file = "neighbor_joining_tree.png") {
  suppressMessages({

    # ️ **Check Input Type & Extract Reference Sequences**
    if (inherits(obj, c("phyloseq", "TreeSummarizedExperiment"))) {
      cat(" Detected object type:", class(obj)[1], "\n Extracting reference sequences...\n")
      ref_sequences <- get_reference_seq(obj)
      if (is.null(ref_sequences)) stop(" Error: No reference sequences found in the object.")
    } else if (is.character(obj) && file.exists(obj)) {
      cat(" Detected FASTA file:", obj, "\n Reading DNA sequences...\n")
      ref_sequences <- Biostrings::readDNAStringSet(obj)
    } else {
      stop(" Invalid input: Must be a phyloseq object, a TreeSummarizedExperiment object, or a FASTA file path.")
    }

    #  **Perform Multiple Sequence Alignment (MSA)**
    cat(" Performing multiple sequence alignment...\n")
    alignment <- msa::msa(ref_sequences)

    # Convert aligned sequences to DNAStringSet
    aligned_sequences <- as(alignment, "DNAStringSet")

    # Convert to phyDat format for distance computation
    phyDat_alignment <- phangorn::phyDat(as(aligned_sequences, "matrix"), type = "DNA")

    # **Compute Distance Matrix using Maximum Likelihood**
    distance_matrix <- phangorn::dist.ml(phyDat_alignment)

    #  **Construct Neighbor-Joining Tree**
    cat(" Constructing Neighbor-Joining tree...\n")
    phylo_tree <- ape::nj(distance_matrix)

    #  **Generate Bootstrap Values**
    cat(" Computing bootstrap support values...\n")
    set.seed(123)  # Ensures reproducibility
    bootstrap_values <- ape::boot.phylo(phylo_tree, as.matrix(aligned_sequences),
                                        FUN = function(x) ape::nj(phangorn::dist.ml(phangorn::phyDat(x, type = "DNA"))),
                                        B = 100)

    # Normalize bootstrap values to percentages
    bootstrap_percentages <- bootstrap_values / 100 * 100

    #  **Polished Plot Aesthetics**
    cat(" Generating plot...\n")
    png(output_file, width = 1000, height = 1000, res = 150)

    par(mar = c(5, 4, 4, 8))  # Increase right margin for labels
    ape::plot.phylo(ape::ladderize(phylo_tree),
                    main = "Neighbor Joining Tree with Bootstrap Values",
                    cex = 1,  # Increase text size
                    edge.width = 1.2,
                    no.margin = TRUE,
                    tip.color = "dodgerblue4")

    ape::nodelabels(round(bootstrap_percentages, 1),
                    cex = 0.9,
                    frame = "none",
                    col = "red4",
                    adj = c(0.5, -0.1))

    dev.off()  # Save plot

    #  **Display Final Tree in Console**
    ape::plot.phylo(ape::ladderize(phylo_tree),
                    main = "Neighbor Joining Tree with Bootstrap Values",
                    cex = 1,
                    edge.width = 1.2,
                    no.margin = TRUE,
                    tip.color = "dodgerblue4")

    ape::nodelabels(round(bootstrap_percentages, 1),
                    cex = 0.9,
                    frame = "none",
                    col = "red4",
                    adj = c(0.5, -0.1))

    cat("\n **Tree successfully saved as:**", output_file, "\n")
  })
}




# Example usage:
# Plot Neighbor-Joining tree with bootstrap values
# combined.fasta <- load_fasta("combined.fasta")
# plot_tree_nj("combind.fasta", output_file = "nj_tree_fasta.png")
