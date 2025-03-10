#' @title Detect Common ASVs and Taxa from Multiple Phyloseq or TSE Objects
#' @description This function identifies the common Amplicon Sequence Variants (ASVs) and
#' common taxa detected across multiple `phyloseq` or `TreeSummarizedExperiment` (TSE) objects.
#' It extracts ASVs and taxa using the package's accessor functions, finds the common ones,
#' and returns a pruned object containing only the shared features.
#'
#' @details Optionally, the results can be saved to CSV and RDS files.
#'
#' @param obj_list A list of `phyloseq` or `TreeSummarizedExperiment` objects.
#' All objects in the list must be of the same class.
#' @param output_common_csv A character string specifying the path to save the pruned ASV/Taxa
#' object as a CSV file. Default is `"common_asvs_taxa.csv"`. Set to `NULL` to disable saving.
#' @param output_common_rds A character string specifying the path to save the pruned ASV/Taxa
#' object as an RDS file. Default is `"common_asvs_taxa.rds"`. Set to `NULL` to disable saving.
#' @param return_as_df A logical indicating whether to return the results as a data frame (`TRUE`)
#' or as the original `phyloseq`/`TSE` object (`FALSE`). Default is `FALSE`.
#'
#' @return A pruned `phyloseq` or `TreeSummarizedExperiment` object containing only common ASVs and taxa.
#' If `return_as_df = TRUE`, a data frame with common ASVs/taxa is returned instead.
#'
#' @importFrom utils write.csv
#' @importFrom phyloseq prune_taxa
#' @importFrom SummarizedExperiment assay
#' @examples
#' \dontrun{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Example with phyloseq objects
#'   common_physeq <- detect_common_asvs_taxa(list(physeq1, physeq2, physeq3))
#'
#'   # Example with TreeSummarizedExperiment objects
#'   common_tse <- detect_common_asvs_taxa(list(tse1, tse2, tse3))
#' }
#' }
#' @export
detect_common_asvs_taxa <- function(obj_list,
                                    output_common_csv = "common_asvs_taxa.csv",
                                    output_common_rds = "common_asvs_taxa.rds",
                                    return_as_df = FALSE) {
  suppressMessages({
    # Validate input
    if (!is.list(obj_list) || length(obj_list) < 2) {
      stop("\U0000274C	obj_list must be a list of at least two phyloseq or TreeSummarizedExperiment objects.")
    }

    # Ensure all objects are of the same type
    obj_classes <- unique(sapply(obj_list, class))
    if (length(obj_classes) > 1) {
      stop("\U0000274C	All objects in obj_list must be of the same class (either all phyloseq or all TreeSummarizedExperiment).")
    }

    # Extract ASVs and taxa using existing accessor functions
    asv_lists <- lapply(obj_list, function(x) tryCatch(rownames(get_otu_table(x)), error = function(e) NULL))
    taxa_lists <- lapply(obj_list, function(x) tryCatch(rownames(get_tax_table(x)), error = function(e) NULL))

    # Remove NULL elements (failed extractions)
    asv_lists <- asv_lists[!sapply(asv_lists, is.null)]
    taxa_lists <- taxa_lists[!sapply(taxa_lists, is.null)]

    # Find common ASVs and taxa
    common_asvs <- Reduce(intersect, asv_lists)
    common_taxa <- Reduce(intersect, taxa_lists)
    common_all <- intersect(common_asvs, common_taxa)  # Ensure ASVs & Taxa are both present

    if (length(common_all) == 0) {
      message("No common ASVs and taxa found.")
      return(NULL)
    }

    # Prune the first obj in the list to retain only common ASVs/taxa
    pruned_obj <- prune_common_taxa(obj_list[[1]], common_all)
    pruned_df <- format_common_data(pruned_obj)

    # Save results if file paths are provided
    if (!is.null(output_common_csv)) utils::write.csv(pruned_df, output_common_csv, row.names = FALSE)
    if (!is.null(output_common_rds)) saveRDS(pruned_obj, file = output_common_rds)

    return(if (return_as_df) pruned_df else pruned_obj)
  })
}

# Example usage:
# results_phy <- detect_common_asvs_taxa(list(physeq1, physeq2))
# results_tse <- detect_common_asvs_taxa(list(tse1, tse2,tse3,tse4))
# detect_common_asvs_taxa(list(tse1, tse2,tse3), output_common_csv
# output.csv", output_common_rds = "output.rds")
