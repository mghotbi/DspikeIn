#' @title Select Important ASVs/OTUs Using Random Forest
#'
#' @description
#' This function selects the most important Amplicon Sequence Variants (ASVs) or Operational Taxonomic Units (OTUs)
#' based on a Random Forest model. If a `TreeSummarizedExperiment (TSE)` is provided, it is first converted to `phyloseq`.
#' The function allows filtering and pruning of taxa before selecting the most important features.
#' Optionally, the selected ASVs/OTUs can be saved as a CSV file.
#'
#' @param physeq A `phyloseq` or `TreeSummarizedExperiment (TSE)` object containing microbiome data.
#' @param response_var A character string specifying the response variable from the sample metadata.
#' @param minlib A numeric value specifying the minimum library size for filtering low-abundance taxa. Default is `15000`.
#' @param prunescale A numeric value specifying the relative abundance threshold for pruning rare OTUs. Default is `0.0001`.
#' @param ntree An integer specifying the number of trees to grow in the Random Forest model. Default is `100`.
#' @param n_top_predictors An integer specifying the number of top ASVs/OTUs to select based on feature importance. Default is `50`.
#' @param output_csv An optional character string specifying the output CSV file name. If `NULL`, no file is saved. Default is `NULL`.
#' @param na_vars A character vector specifying metadata variables to check for missing values (`NA`). If `NULL`, only `response_var` is checked.
#'
#' @return
#' Returns a pruned `phyloseq` or `TreeSummarizedExperiment (TSE)` object containing only the selected ASVs/OTUs.
#' If the input is `TSE`, the output is converted back to `TSE`.
#'
#' @importFrom phyloseq subset_taxa sample_data tax_table sample_sums prune_samples prune_taxa otu_table psmelt taxa_sums
#' @importFrom dplyr arrange desc
#' @importFrom randomForest randomForest importance
#' @importFrom utils write.csv
#'
#' @examples
#' \donttest{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Perform Random Forest feature selection
#'   rf_physeq <- RandomForest_selected_ASVs(
#'     physeq_16SOTU,
#'     response_var = "Host.genus",
#'     na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet")
#'   )
#'
#'   # Load TreeSummarizedExperiment (TSE) object
#'   tse_16SOTU <- convert_phyloseq_to_tse(physeq_16SOTU)
#'
#'   # Perform Random Forest feature selection on TSE object
#'   rf_tse <- RandomForest_selected_ASVs(
#'     tse_16SOTU,
#'     response_var = "Host.genus",
#'     na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet")
#'   )
#' }
#' }
#'
#' @export
RandomForest_selected_ASVs <- function(physeq, response_var, minlib = 15000, prunescale = 0.0001,
                                       ntree = 100, n_top_predictors = 100, output_csv = NULL,
                                       na_vars = NULL) {
  suppressMessages({

    #  Detect if input is TSE and convert to phyloseq
    is_TSE <- inherits(physeq, "TreeSummarizedExperiment")
    if (is_TSE) {
      physeq <- convert_tse_to_phyloseq(physeq)
    }

    #  Check if response_var exists in sample_data
    if (!response_var %in% colnames(phyloseq::sample_data(physeq))) {
      stop(paste("Response variable", response_var, "not found in sample data."))
    }

    #  Ensure na_vars includes response_var
    if (!is.null(na_vars)) {
      na_vars <- unique(c(response_var, na_vars))
    } else {
      na_vars <- response_var
    }

    #  Remove samples with NA in specified variables
    sample_data_df <- data.frame(phyloseq::sample_data(physeq))
    for (var in na_vars) {
      sample_data_df <- sample_data_df[!is.na(sample_data_df[[var]]), ]
    }
    phyloseq::sample_data(physeq) <- phyloseq::sample_data(sample_data_df)

    #  Convert categorical variables to factors
    for (var in na_vars) {
      phyloseq::sample_data(physeq)[[var]] <- as.factor(phyloseq::sample_data(physeq)[[var]])
    }

    #  Ensure response variable has at least two levels
    if (length(unique(sample_data_df[[response_var]])) < 2) {
      stop("The response variable needs at least two classes for classification.")
    }

    #  Remove taxa with any NA values in tax_table
    physeq <- phyloseq::subset_taxa(physeq, apply(phyloseq::tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))

    #  Remove samples with zero reads after subsetting
    physeq <- phyloseq::prune_samples(phyloseq::sample_sums(physeq) > 0, physeq)

    #  Prune out rare OTUs using prunescale
    tax.mean <- phyloseq::taxa_sums(physeq) / phyloseq::nsamples(physeq)
    s.prune <- phyloseq::prune_taxa(tax.mean > prunescale * minlib, physeq)

    #  Replace empty/NA values in tax_table with "unidentified"
    phyloseq::tax_table(s.prune) <- replace(phyloseq::tax_table(s.prune), is.na(phyloseq::tax_table(s.prune)), "unidentified")
    phyloseq::tax_table(s.prune) <- replace(phyloseq::tax_table(s.prune), phyloseq::tax_table(s.prune) == "", "unidentified")

    predictors <- t(phyloseq::otu_table(s.prune))

    #  Prepare response variable
    response <- as.factor(phyloseq::sample_data(s.prune)[[response_var]])

    #  Combine into one data frame
    rf.data <- data.frame(response, predictors)

    #  emove rows with NA values from rf.data
    rf.data <- rf.data[complete.cases(rf.data), ]

    #  Check again if response variable has at least two classes
    if (length(unique(rf.data$response)) < 2) {
      stop("\U0000274C The response variable needs at least two classes for classification after removing NAs.")
    }

    #  Train Random Forest model
    set.seed(2)
    sal.classify <- randomForest::randomForest(response ~ ., data = rf.data, ntree = ntree, keep.forest = FALSE, proximity = FALSE, importance = TRUE)

    #  Extract variable importance
    imp <- randomForest::importance(sal.classify)
    imp_df <- data.frame(predictors = rownames(imp), imp)
    imp_df <- dplyr::arrange(imp_df, dplyr::desc(MeanDecreaseGini))
    imp_df$predictors <- factor(imp_df$predictors, levels = imp_df$predictors)

    #  Select the top predictors
    imp_top <- imp_df[1:n_top_predictors, ]

    # Subset selected ASVs/OTUs from phyloseq object
    otunames <- imp_top$predictors
    r <- rownames(phyloseq::tax_table(s.prune)) %in% otunames
    selected_phyloseq <- phyloseq::prune_taxa(r, s.prune)

    if (sum(r) > 0) {
      # Save selected ASVs to CSV **only if output_csv is provided**
      if (!is.null(output_csv)) {
        pm <- phyloseq::psmelt(selected_phyloseq)
        utils::write.csv(pm, file = output_csv, row.names = FALSE)
        cat("\U0001F4BE Selected ASVs saved to:", output_csv, "\n")
      }

      # Convert back to TSE if needed
      if (is_TSE) {
        selected_phyloseq <- convert_phyloseq_to_tse(selected_phyloseq)
      }

      return(selected_phyloseq)
    } else {
      print("No ASVs selected.")

      # Convert back to TSE if needed (even if no ASVs selected)
      if (is_TSE) {
        return(convert_phyloseq_to_tse(selected_phyloseq))
      }

      return(NULL)
    }
  })
}

# Example usage:
# rf_physeq <- RandomForest_selected_ASVs(physeq_16SOTU, response_var = "Host.genus",
# na_vars = c("Habitat", "Ecoregion.III", "Host.genus", "Diet"))
# saveRDS(rf_physeq,"rf_physeq.rds")
# rf_physeq@tax_table
