#' @title Validate Spike-In Clade Consistency with NJ Tree and Bootstrap
#'
#' @description
#' Validates whether sample spike-in sequences form a monophyletic clade
#' with known reference spike-in(s) using a Neighbor-Joining (NJ) tree
#' with Jukes-Cantor (JC69) correction and bootstrap support.
#'
#' Sequences are aligned with \code{DECIPHER::AlignSeqs()}, converted to a
#' \code{phyDat} object, and a JC69 distance matrix is used to build an NJ
#' tree. The tree is rooted on the first reference sequence, and the
#' bootstrap support reported is the support of the node joining all sample
#' sequences (their most recent common ancestor), not the maximum support
#' found anywhere in the tree.
#'
#' This function produces:
#' \itemize{
#'   \item A bootstrap-annotated NJ tree
#'   \item A boxplot comparing terminal branch lengths
#'   \item A histogram of patristic distances to the reference(s)
#' }
#' If \code{output_prefix} is provided, outputs are saved to disk. In all
#' cases the plots are also returned as recorded plots.
#'
#' @details
#' With a single reference sequence, monophyly of the sample sequences is
#' guaranteed by construction (the tree is rooted on that reference), so
#' \code{monophyly} is only informative, and \code{clade_bootstrap} is only
#' computed, when two or more reference sequences are supplied. The
#' patristic distances to the reference are informative in either case.
#'
#' @param reference_fasta Character. Path to FASTA file of reference spike-ins.
#' @param sample_fasta Character. Path to FASTA file of sample spike-ins.
#' @param bootstrap Integer. Number of bootstrap replicates (default = 100).
#' @param output_prefix Character or NULL. File prefix for saving output.
#'
#' @return An invisible list with:
#' \describe{
#'   \item{tree}{NJ tree rooted on the first reference (class \code{phylo})}
#'   \item{monophyly}{TRUE if sample spike-ins form a clade}
#'   \item{clade_bootstrap}{Bootstrap support (\%) of the sample clade,
#'     or \code{NA} if the samples are not monophyletic or fewer than two
#'     reference or two sample sequences are supplied}
#'   \item{node_support}{Bootstrap support (\%) for every internal node
#'     (\code{NA} for the root)}
#'   \item{branch_stats}{Terminal branch length summary by group}
#'   \item{branch_table}{Terminal branch length for every tip, taken from
#'     the unrooted NJ tree}
#'   \item{patristic_distances}{Patristic distance matrix}
#'   \item{tree_plot}{Recorded tree plot}
#'   \item{branch_boxplot}{Recorded boxplot}
#'   \item{patristic_histogram}{Recorded histogram}
#'   \item{summary_text}{Text summary}
#'   \item{alignment}{Aligned sequences (\code{DNAStringSet})}
#'   \item{aln_phydat}{Alignment converted to \code{phyDat}}
#'   \item{distance_matrix}{JC69 distance matrix}
#' }
#'
#' @examples
#' ref_fasta <- system.file("extdata", "Ref.fasta", package = "DspikeIn")
#' sample_fasta <- system.file("extdata", "Sample.fasta", package = "DspikeIn")
#' result <- validate_spikein_clade(ref_fasta, sample_fasta, bootstrap = 20)
#' result$summary_text
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom DECIPHER AlignSeqs
#' @importFrom phangorn phyDat dist.ml bootstrap.phyDat
#' @importFrom ape nj root ladderize is.monophyletic getMRCA cophenetic.phylo
#'   prop.clades plot.phylo nodelabels edgelabels tiplabels write.tree
#' @importFrom grDevices png pdf dev.off dev.control recordPlot
#' @importFrom graphics legend boxplot hist
#' @importFrom stats aggregate sd
#' @importFrom utils write.csv
#' @export
validate_spikein_clade <- function(reference_fasta,
                                   sample_fasta,
                                   bootstrap = 100,
                                   output_prefix = NULL) {
  ## ---- Input checks ------------------------------------------------------
  if (!is.character(reference_fasta) || length(reference_fasta) != 1L ||
    !file.exists(reference_fasta)) {
    stop("'reference_fasta' must be the path to an existing FASTA file.")
  }
  if (!is.character(sample_fasta) || length(sample_fasta) != 1L ||
    !file.exists(sample_fasta)) {
    stop("'sample_fasta' must be the path to an existing FASTA file.")
  }
  if (!is.numeric(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap) ||
    bootstrap < 1) {
    stop("'bootstrap' must be a single positive integer.")
  }
  bootstrap <- as.integer(bootstrap)
  if (!is.null(output_prefix) &&
    (!is.character(output_prefix) || length(output_prefix) != 1L)) {
    stop("'output_prefix' must be NULL or a single character string.")
  }

  ## ---- Read sequences ----------------------------------------------------
  message("Loading reference sequences...")
  ref_seqs <- Biostrings::readDNAStringSet(reference_fasta)
  message("Loading sample sequences...")
  sample_seqs <- Biostrings::readDNAStringSet(sample_fasta)
  if (length(ref_seqs) < 1L) stop("Reference FASTA is empty.")
  if (length(sample_seqs) < 1L) stop("Sample FASTA is empty.")

  # Clean headers once, up front: use the 'sample=' tag when present and
  # make every label unique (phyDat / dist.ml require unique names).
  clean_label <- function(header) {
    ifelse(grepl("sample=", header, fixed = TRUE),
      sub(".*sample=([^;[:space:]]+).*", "\\1", header),
      header
    )
  }
  ref_labels_raw <- clean_label(names(ref_seqs))
  sample_labels_raw <- clean_label(names(sample_seqs))
  all_labels <- make.unique(c(ref_labels_raw, sample_labels_raw), sep = "_")
  ref_labels <- all_labels[seq_along(ref_seqs)]
  sample_labels <- all_labels[-seq_along(ref_seqs)]

  combined_seqs <- c(ref_seqs, sample_seqs)
  names(combined_seqs) <- all_labels
  if (length(combined_seqs) < 3L) {
    stop("At least three sequences in total are required to build a tree.")
  }
  message("Total sequences combined: ", length(combined_seqs))

  ## ---- Alignment ---------------------------------------------------------
  message("Performing multiple sequence alignment (DECIPHER::AlignSeqs)...")
  alignment <- DECIPHER::AlignSeqs(combined_seqs,
    processors = 1L, verbose = FALSE
  )
  aln_chars <- strsplit(as.character(alignment), "", fixed = TRUE)
  aln_mat <- do.call(rbind, aln_chars)
  rownames(aln_mat) <- names(alignment)
  aln_phydat <- phangorn::phyDat(aln_mat, type = "DNA")

  ## ---- Distance + NJ tree ------------------------------------------------
  message("Computing JC69 distance matrix...")
  dist_jc69 <- phangorn::dist.ml(aln_phydat, model = "JC69")
  nj_tree <- ape::nj(dist_jc69)
  tree <- ape::root(nj_tree, outgroup = ref_labels[1], resolve.root = TRUE)
  tree <- ape::ladderize(tree)
  n_tip <- length(tree$tip.label)

  is_ref <- tree$tip.label %in% ref_labels
  is_clade <- if (length(sample_labels) >= 2L) {
    ape::is.monophyletic(tree, sample_labels)
  } else {
    TRUE
  }

  ## ---- Bootstrap ---------------------------------------------------------
  message("Performing bootstrap (n = ", bootstrap, ") ...")
  bs_trees <- phangorn::bootstrap.phyDat(
    aln_phydat,
    FUN = function(x) ape::nj(phangorn::dist.ml(x, model = "JC69")),
    bs = bootstrap
  )
  # Bipartition-based support, since the NJ replicates are unrooted.
  node_support <- ape::prop.clades(tree, bs_trees, rooted = FALSE) /
    bootstrap * 100

  # Clade support is only testable with >= 2 references and >= 2 samples;
  # with a single reference the sample clade is guaranteed by the rooting.
  clade_testable <- length(ref_labels) >= 2L && length(sample_labels) >= 2L
  clade_bootstrap <- NA_real_
  if (is_clade && clade_testable) {
    mrca <- ape::getMRCA(tree, sample_labels)
    clade_bootstrap <- round(node_support[mrca - n_tip], 1)
  }

  ## ---- Terminal branch lengths (in tip order) ----------------------------
  # Taken from the unrooted NJ tree: rooting on a reference splits its
  # terminal branch and would artificially shorten it.
  nj_tip_edge <- match(
    match(tree$tip.label, nj_tree$tip.label),
    nj_tree$edge[, 2]
  )
  branch_table <- data.frame(
    tip = tree$tip.label,
    group = ifelse(is_ref, "Reference", "Sample"),
    branch_length = nj_tree$edge.length[nj_tip_edge],
    stringsAsFactors = FALSE
  )
  stat_summary <- stats::aggregate(branch_length ~ group,
    data = branch_table,
    FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
  )

  patristic_dist <- ape::cophenetic.phylo(tree)
  dist_to_ref <- patristic_dist[sample_labels, ref_labels, drop = FALSE]
  mean_dist_to_ref <- rowMeans(dist_to_ref)

  ## ---- Plots -------------------------------------------------------------
  group_cols <- c(Reference = "#3F37C9", Sample = "#FF5722")

  draw_enhanced_tree <- function() {
    ape::plot.phylo(tree,
      type = "phylogram", edge.width = 1.5, cex = 0.7,
      tip.color = group_cols[branch_table$group],
      main = "Spike-In Clade Validation Tree (NJ, JC69)"
    )
    legend("topright",
      legend = names(group_cols), col = group_cols,
      pch = 19, bty = "n"
    )
    ape::nodelabels(
      text = ifelse(is.na(node_support), "", round(node_support)),
      frame = "none",
      adj = c(1.1, -0.3), cex = 0.6
    )
    ape::edgelabels(
      text = round(tree$edge.length, 3), frame = "none",
      adj = c(0.5, 1.3), cex = 0.5, col = "gray40"
    )
    tip_d <- rep("", n_tip)
    tip_d[match(sample_labels, tree$tip.label)] <-
      paste0("d=", round(mean_dist_to_ref, 3))
    ape::tiplabels(tip_d, frame = "none", adj = c(-0.1, 1.5), cex = 0.5)
  }

  draw_boxplot <- function() {
    boxplot(branch_length ~ group,
      data = branch_table,
      col = group_cols[sort(unique(branch_table$group))],
      main = "Terminal Branch Length Comparison", ylab = "Branch Length"
    )
  }

  draw_hist <- function() {
    hist(as.vector(dist_to_ref),
      main = "Patristic Distances to Reference",
      xlab = "Patristic Distance", col = "gray80", border = "white"
    )
  }

  # Record a plot on an off-screen device so the function never opens
  # windows or leaves stray Rplots.pdf files behind.
  record_plot <- function(plot_fn) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::dev.control(displaylist = "enable")
    plot_fn()
    grDevices::recordPlot()
  }

  tree_plot <- record_plot(draw_enhanced_tree)
  branch_boxplot <- record_plot(draw_boxplot)
  patristic_histogram <- record_plot(draw_hist)

  if (!is.null(output_prefix)) {
    save_plot <- function(open_dev, plot_fn) {
      open_dev()
      on.exit(grDevices::dev.off(), add = TRUE)
      plot_fn()
    }
    save_plot(function() {
      grDevices::png(paste0(output_prefix, "_tree.png"),
        width = 2000, height = 2000, res = 300
      )
    }, draw_enhanced_tree)
    save_plot(function() {
      grDevices::pdf(paste0(output_prefix, "_branch_lengths.pdf"))
    }, draw_boxplot)
    save_plot(function() {
      grDevices::pdf(paste0(output_prefix, "_patristic_distances_hist.pdf"))
    }, draw_hist)

    ape::write.tree(tree, file = paste0(output_prefix, ".nwk"))
    utils::write.csv(as.matrix(patristic_dist),
      file = paste0(output_prefix, "_patristic_distances.csv")
    )
  }

  ## ---- Summary -----------------------------------------------------------
  dist_range <- round(range(dist_to_ref), 3)
  dist_mean <- round(mean(dist_to_ref), 3)
  clade_part <- if (!clade_testable) {
    paste0(
      length(sample_labels), " sample sequence(s) compared with ",
      length(ref_labels), " reference(s); clade support is not testable ",
      "with fewer than two references or two samples. "
    )
  } else {
    paste0(
      length(sample_labels), " sample sequences ",
      if (is_clade) "formed" else "did NOT form",
      " a single clade relative to ", length(ref_labels),
      " references (bootstrap = ",
      if (is.na(clade_bootstrap)) "NA" else paste0(clade_bootstrap, "%"),
      "). "
    )
  }
  summary_text <- paste0(
    "Validation Result: ", clade_part,
    "The mean patristic distance to the reference is ", dist_mean,
    ", range [", dist_range[1], ", ", dist_range[2], "]."
  )
  message(summary_text)

  invisible(list(
    tree = tree,
    monophyly = is_clade,
    clade_bootstrap = clade_bootstrap,
    node_support = node_support,
    branch_stats = stat_summary,
    branch_table = branch_table,
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
