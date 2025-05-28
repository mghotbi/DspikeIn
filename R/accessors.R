#' @title Extract OTU Tax Metadata from Object
#' @description Retrieves the OTU table from a `phyloseq` or `TreeSummarizedExperiment` object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A matrix containing OTU count data.
#' @importFrom phyloseq otu_table
#' @importFrom SummarizedExperiment assay
#' @details This is a thin wrapper that unifies access to `assay()` (for TreeSummarizedExperiment) and `otu_table()` (for phyloseq) into a single function for format-agnostic use.
#' @seealso \code{\link[phyloseq]{otu_table}}, \code{\link[SummarizedExperiment]{assay}}
#' @export
get_otu_table <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    return(as(otu_table(obj), "matrix"))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(SummarizedExperiment::assay(obj))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}

#' @title Extract Taxonomy Table
#' @description Retrieves the taxonomy table from an object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A data frame containing taxonomy annotations.
#' @importFrom phyloseq tax_table
#' @importFrom SummarizedExperiment rowData
#' @details This is a thin wrapper around `rowData()` or `tax_table()` to ensure compatibility across data structures.
#' @seealso \code{\link[phyloseq]{tax_table}}, \code{\link[SummarizedExperiment]{rowData}}

#' @export
get_tax_table <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    return(as.data.frame(tax_table(obj)))
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::rowData(obj)))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}

#' @title Extract Sample Data
#' @description Retrieves the sample metadata.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A data frame with sample metadata.
#' @importFrom phyloseq sample_data
#' @importFrom microbiome meta
#' @importFrom SummarizedExperiment colData
#' @details Thin wrapper over `colData()` or `sample_data()` for unified metadata access across formats.
#' @seealso \code{\link[phyloseq]{sample_data}}, \code{\link[SummarizedExperiment]{colData}}, \code{\link[microbiome]{meta}}
#' @export
get_sample_data <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    return(as.data.frame(microbiome::meta(obj))) # Preferred method for phyloseq
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(obj)))
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}

#' @title Extract Phylogenetic Tree
#' @description Retrieves the phylogenetic tree.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A phylogenetic tree object.
#' @importFrom phyloseq phy_tree
#' @importFrom TreeSummarizedExperiment rowTree
#' @details This function wraps `phy_tree()` and `rowTree()` for compatibility with both frameworks.
#' @seealso \code{\link[phyloseq]{phy_tree}}, \code{\link[TreeSummarizedExperiment]{rowTree}}
#' @export
get_phy_tree <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    tree <- phyloseq::phy_tree(obj)
    if (is.null(tree)) stop("No phylogenetic tree found in phyloseq object.")
    return(tree)
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    if (!is.null(TreeSummarizedExperiment::rowTree(obj))) {
      return(TreeSummarizedExperiment::rowTree(obj))
    } else {
      stop("No tree found in rowTree(obj).")
    }
  } else {
    stop("Unsupported object type. Must be phyloseq or TreeSummarizedExperiment.")
  }
}


#' @title Extract Reference Sequences
#' @description Retrieves the reference sequences (if available) from an object.
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A `DNAStringSet` object containing the reference sequences.
#' @importFrom phyloseq refseq
#' @importFrom S4Vectors metadata
#' @details Thin wrapper around `refseq()` or `metadata(obj)$refseq` for unified sequence retrieval.
#' @seealso \code{\link[phyloseq]{refseq}}, \code{\link[S4Vectors]{metadata}}
#' @export
get_reference_seq <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    if (!is.null(phyloseq::refseq(obj))) {
      return(phyloseq::refseq(obj))
    } else {
      stop(" No reference sequences found in the phyloseq object.")
    }
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    if ("refseq" %in% names(S4Vectors::metadata(obj))) {
      return(S4Vectors::metadata(obj)$refseq)
    } else {
      stop(" No reference sequences found in the TreeSummarizedExperiment object.")
    }
  } else {
    stop("Unsupported object type: must be phyloseq or TreeSummarizedExperiment.")
  }
}
