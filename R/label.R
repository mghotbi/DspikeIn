#' @title Label Taxonomic Ranks by Hashcode
#'
#' @description
#' Labels ASVs/OTUs in a `phyloseq` object using a named vector mapping hashcodes
#' to known taxonomy labels (e.g., spike-in species, genera, or families). This is
#' especially useful for clearly labeling synthetic controls or mock taxa.
#'
#' If the specified taxonomic rank (e.g., "Genus" or "Family") does not exist in the
#' `tax_table`, it will be added and filled with `NA`, followed by inserting your custom labels.
#'
#' @param obj A [`phyloseq::phyloseq`] object with a populated `tax_table()`.
#' @param hashcode_label_map A named character vector where names are ASV hashcodes
#'   (i.e., row names of the tax table), and values are the taxonomy labels to assign.
#' @param tax_rank A character string specifying the taxonomic rank to label
#'   (e.g., `"Species"`, `"Genus"`, `"Family"`). Default is `"Species"`.
#'
#' @return A modified `phyloseq` object with updated taxonomic labels at the specified rank.
#'
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   library(phyloseq)
#'
#'   # Create dummy tax_table with hashcodes
#'   tax_mat <- matrix(
#'     data = c(
#'       "Bacteria", NA,
#'       "Bacteria", NA,
#'       "Bacteria", NA,
#'       "Bacteria", NA
#'     ),
#'     nrow = 4,
#'     dimnames = list(
#'       c(
#'         "8ac7ad6e4b6501eb143d97f10bcc2b6d",
#'         "5a92565231c6df8f58871c0df2d1a12a",
#'         "8bf5a7b04cb725bc3d2627b971eb03fb",
#'         "7ae171e44f46ddadcbb53ee6b34a483b"
#'       ),
#'       c("Kingdom", "Species")
#'     )
#'   )
#'
#'   ps <- phyloseq::phyloseq(phyloseq::tax_table(tax_mat))
#'
#'   # Define hash-to-label map
#'   hash_map <- c(
#'     "8ac7ad6e4b6501eb143d97f10bcc2b6d" = "MockGenus_A",
#'     "5a92565231c6df8f58871c0df2d1a12a" = "MockGenus_B"
#'   )
#'
#'   # Label genus (automatically adds "Genus" column if missing)
#'   ps <- label(ps, hash_map, tax_rank = "Genus")
#'
#'   # Check updated taxonomy
#'   phyloseq::tax_table(ps)[, "Genus"]
#' }
#'
#' @seealso [phyloseq::tax_table()]
#' @importFrom phyloseq tax_table<- tax_table
#' @export
label <- function(obj, hashcode_label_map, tax_rank = "Species") {
  tax_df <- as.data.frame(phyloseq::tax_table(obj))

  # Add the taxonomic rank if missing
  if (!tax_rank %in% colnames(tax_df)) {
    warning("Taxonomic rank '", tax_rank, "' not found in tax_table. Adding it and filling with NA.")
    tax_df[[tax_rank]] <- NA
  }

  # Check that hashcodes exist
  missing <- setdiff(names(hashcode_label_map), rownames(tax_df))
  if (length(missing) > 0) {
    stop(
      "The following hashcodes were not found in the tax_table: ",
      paste(missing, collapse = ", ")
    )
  }

  # Update tax_table entries at the specified rank
  tax_df[names(hashcode_label_map), tax_rank] <- hashcode_label_map

  # Update the tax_table in the phyloseq object
  phyloseq::tax_table(obj) <- phyloseq::tax_table(as.matrix(tax_df))

  message("Updated taxonomic rank '", tax_rank, "' for matched ASV hashcodes.")
  return(obj)
}
