#' @title Convert a Phyloseq or TSE Object into a Long-Format Data Frame
#' @description Converts a `phyloseq` or `TreeSummarizedExperiment` object into
#' a long-format data frame, similar to `phyloseq::psmelt()`, for compatibility
#' with visualization functions such as `alluvial_plot()`.
#'
#' @param obj A `phyloseq` or `TreeSummarizedExperiment` object.
#' @return A long-format `data.frame` containing taxonomic, abundance, and sample metadata.
#'
#' @importFrom phyloseq psmelt tax_table otu_table sample_data
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate select rename bind_cols
#' @importFrom tidyr pivot_longer
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Turn it to the long format
#'   melted_ph <- get_long_format_data(physeq_16SOTU)
#'   melted <- get_long_format_data(tse_16SOTU)
#' }
#' @export
get_long_format_data <- function(obj) {
  if (inherits(obj, "phyloseq")) {
    message("Converting phyloseq object to long format using psmelt()...")
    return(phyloseq::psmelt(obj)) # Use psmelt directly for phyloseq
  } else if (inherits(obj, "TreeSummarizedExperiment")) {
    message("Converting TreeSummarizedExperiment object to long format...")

    # Extract abundance data
    abundance_data <- as.data.frame(SummarizedExperiment::assay(obj))
    abundance_data <- tibble::rownames_to_column(abundance_data, "ASV")

    # Extract taxonomic annotations
    tax_data <- as.data.frame(SummarizedExperiment::rowData(obj))
    tax_data <- tibble::rownames_to_column(tax_data, "ASV")

    # Extract sample metadata
    sample_data <- as.data.frame(SummarizedExperiment::colData(obj))
    sample_data <- tibble::rownames_to_column(sample_data, "Sample")

    # Reshape abundance data to long format
    long_data <- tidyr::pivot_longer(
      abundance_data,
      cols = -ASV,
      names_to = "Sample",
      values_to = "Abundance"
    )

    # Merge taxonomy and metadata
    long_data <- dplyr::left_join(long_data, tax_data, by = "ASV")
    long_data <- dplyr::left_join(long_data, sample_data, by = "Sample")

    return(long_data)
  } else {
    stop("Unsupported object type. Input must be a phyloseq or TreeSummarizedExperiment object.")
  }
}

# Usage Example
# melted<-get_long_format_data(tse_F)
# melted_ph<-get_long_format_data(phy_M_M)
