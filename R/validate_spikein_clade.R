#' @title Validate Spike-In Clade Consistency with NJ Tree and Bootstrap
#'
#' @description
#' Validates whether sample spike-in sequences form a monophyletic clade
#' with known reference spike-in(s) using a Neighbor-Joining (NJ) tree
#' with Jukes-Cantor correction and bootstrap support.
#'
#' This function produces:
#' \itemize{
#'   \item A bootstrap-annotated NJ tree
#'   \item A boxplot comparing branch lengths
#'   \item A histogram of patristic distances
#' }
#' If \code{output_prefix} is provided, outputs are saved. Otherwise, they
#' are shown interactively and recorded in-memory.
#'
#' @param reference_fasta Character. Path to FASTA file of reference spike-ins.
#' @param sample_fasta Character. Path to FASTA file of sample spike-ins.
#' @param bootstrap Integer. Number of bootstrap replicates (default = 100).
#' @param output_prefix Character or NULL. File prefix for saving output.
#'
#' @return A list with:
#' \describe{
#'   \item{tree}{NJ tree (class \code{phylo})}
#'   \item{monophyly}{TRUE if sample spike-ins form a clade}
#'   \item{clade_bootstrap}{Bootstrap support percentage}
#'   \item{branch_stats}{Branch length summary}
#'   \item{patristic_distances}{Patristic distance matrix}
#'   \item{tree_plot}{Tree plot object}
#'   \item{branch_boxplot}{Boxplot object}
#'   \item{patristic_histogram}{Histogram object}
#'   \item{summary_text}{Text summary}
#'   \item{alignment}{Multiple sequence alignment (MSA)}
#'   \item{aln_phydat}{Alignment converted to \code{phyDat}}
#'   \item{distance_matrix}{JC69 distance matrix}
#' }
#'
#' @examples
#' ref_fasta <- system.file("extdata", "Ref.fasta", package = "DspikeIn")
#' sample_fasta <- system.file("extdata", "Sample.fasta", package = "DspikeIn")
#' result <- validate_spikein_clade(ref_fasta, sample_fasta)
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom msa msa msaConvert
#' @importFrom phangorn phyDat dist.ml bootstrap.phyDat
#' @importFrom ape nj ladderize is.monophyletic cophenetic.phylo prop.clades plot.phylo nodelabels tiplabels write.tree
#' @importFrom grDevices png pdf dev.off dev.new
#' @importFrom graphics legend boxplot hist
#' @importFrom utils write.csv
#' @export
validate_spikein_clade <- function(reference_fasta,
                                   sample_fasta,
                                   bootstrap = 100,
                                   output_prefix = NULL) {
  message("Loading reference sequences...")
  ref_seqs <- Biostrings::readDNAStringSet(reference_fasta)

  message("Loading sample sequences...")
  sample_seqs <- Biostrings::readDNAStringSet(sample_fasta)
  if (length(ref_seqs) < 1) stop("Reference FASTA is empty.")
  if (length(sample_seqs) < 1) stop("Sample FASTA is empty.")

  combined_seqs <- c(ref_seqs, sample_seqs)
  message("Total sequences combined: ", length(combined_seqs))

  if (!requireNamespace("msa", quietly = TRUE)) stop("Please install the 'msa' package.")

  message("Performing multiple sequence alignment...")
  alignment <- msa::msa(combined_seqs)
  aln_phydat <- msa::msaConvert(alignment, type = "phangorn::phyDat")

  message("Computing JC69 distance matrix...")
  dist_jc69 <- phangorn::dist.ml(aln_phydat)
  tree <- ape::ladderize(ape::nj(dist_jc69))

  # Clean labels and map sample IDs
  tip_map <- vapply(
    names(combined_seqs),
    function(header) {
      if (grepl("sample=", header)) {
        sub(".*sample=([^;\\s]+).*", "\\1", header)
      } else {
        header
      }
    },
    FUN.VALUE = character(1)
  )
  names(tip_map) <- names(combined_seqs)
  tree$tip.label <- tip_map[tree$tip.label]

  ref_tips <- tip_map[names(ref_seqs)]
  is_ref <- tree$tip.label %in% ref_tips
  is_clade <- ape::is.monophyletic(tree, which(!is_ref))

  message("Performing bootstrap (n = ", bootstrap, ") ...")
  bs_tree <- phangorn::bootstrap.phyDat(
    aln_phydat,
    FUN = function(x) ape::nj(phangorn::dist.ml(x)),
    bs = bootstrap
  )

  clade_nodes <- ape::prop.clades(tree, bs_tree) / bootstrap * 100
  clade_bootstrap <- if (is_clade) round(max(clade_nodes, na.rm = TRUE), 1) else NA

  terminal_edges <- which(tree$edge[, 2] <= length(tree$tip.label))
  branch_lengths <- tree$edge.length[terminal_edges]
  group <- ifelse(is_ref, "Reference", "Sample")

  branch_table <- data.frame(
    tip = tree$tip.label,
    group = group,
    branch_length = branch_lengths
  )

  stat_summary <- aggregate(branch_length ~ group,
    data = branch_table,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )

  patristic_dist <- ape::cophenetic.phylo(tree)
  dist_to_ref <- patristic_dist[tree$tip.label[!is_ref], tree$tip.label[is_ref], drop = FALSE]

  # Plot builder
  draw_enhanced_tree <- function() {
    ape::plot.phylo(tree,
      type = "phylogram", edge.width = 1.5, cex = 0.9,
      tip.color = ifelse(is_ref, "#3F37C9", "#FF5722"),
      main = "Spike-In Clade Validation Tree (JC69)"
    )
    legend("topright",
      legend = c("Reference", "Sample"),
      col = c("#3F37C9", "#FF5722"), pch = 19, bty = "n"
    )
    ape::nodelabels(
      text = round(clade_nodes, 1), frame = "none",
      adj = c(1.1, 0.2), cex = 0.7
    )
    ape::nodelabels(
      text = paste0("L=", round(tree$edge.length, 3)),
      frame = "none", adj = c(-0.2, -0.6), cex = 0.6, col = "gray40"
    )
    label_values <- if (ncol(dist_to_ref) == 1) {
      round(dist_to_ref[, 1], 3)
    } else {
      round(apply(dist_to_ref, 1, mean), 3)
    }
    ape::tiplabels(paste0("d=", label_values), frame = "none", adj = c(1, -0.5), cex = 0.7)
  }

  ## === Plot Management ===
  record_plot <- function(plot_fn) {
    grDevices::dev.new()
    plot_fn()
    p <- recordPlot()
    dev.off()
    p
  }

  tree_plot <- branch_boxplot <- patristic_histogram <- NULL

  if (is.null(output_prefix)) {
    tree_plot <- record_plot(draw_enhanced_tree)

    branch_boxplot <- record_plot(function() {
      boxplot(branch_length ~ group,
        data = branch_table,
        col = c("#3F37C9", "#FF5722"),
        main = "Branch Length Comparison", ylab = "Branch Length"
      )
    })

    patristic_histogram <- record_plot(function() {
      hist(dist_to_ref,
        main = "Patristic Distances to Reference",
        xlab = "Patristic Distance"
      )
    })
  } else {
    grDevices::png(paste0(output_prefix, "_tree.png"), width = 2000, height = 2000, res = 500)
    draw_enhanced_tree()
    dev.off()

    grDevices::pdf(paste0(output_prefix, "_branch_lengths.pdf"))
    boxplot(branch_length ~ group,
      data = branch_table,
      col = c("#3F37C9", "#FF5722"),
      main = "Branch Length Comparison", ylab = "Branch Length"
    )
    dev.off()

    grDevices::pdf(paste0(output_prefix, "_patristic_distances_hist.pdf"))
    hist(dist_to_ref,
      main = "Patristic Distances to Reference",
      xlab = "Patristic Distance"
    )
    dev.off()

    # Also record interactive versions
    tree_plot <- record_plot(draw_enhanced_tree)
    branch_boxplot <- record_plot(function() {
      boxplot(branch_length ~ group,
        data = branch_table,
        col = c("#3F37C9", "#FF5722"),
        main = "Branch Length Comparison", ylab = "Branch Length"
      )
    })
    patristic_histogram <- record_plot(function() {
      hist(dist_to_ref,
        main = "Patristic Distances to Reference",
        xlab = "Patristic Distance"
      )
    })

    ape::write.tree(tree, file = paste0(output_prefix, ".nwk"))
    utils::write.csv(as.matrix(patristic_dist),
      file = paste0(output_prefix, "_patristic_distances.csv")
    )
  }

  # Build summary
  ref_count <- sum(is_ref)
  sample_count <- sum(!is_ref)
  dist_range <- round(range(dist_to_ref), 3)
  dist_mean <- round(mean(dist_to_ref), 3)
  summary_text <- paste0(
    "Validation Result: ", sample_count,
    " sample sequences formed a clade with ", ref_count,
    " reference(s) (bootstrap = ", clade_bootstrap, "%). ",
    "The mean patristic distance to the reference is ", dist_mean,
    ", range [", dist_range[1], ", ", dist_range[2], "]."
  )
  message(summary_text)

  invisible(list(
    tree = tree,
    monophyly = is_clade,
    clade_bootstrap = clade_bootstrap,
    branch_stats = stat_summary,
    patristic_distances = patristic_dist,
    tree_plot = tree_plot,
    branch_boxplot = branch_boxplot,
    patristic_histogram = patristic_histogram,
    summary_text = summary_text,
    alignment = alignment,
    aln_phydat = aln_phydat,
    distance_matrix = dist_jc69
  ))
}
