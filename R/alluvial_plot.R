#' @title Generate an Alluvial Plot for Microbiome Data
#' @note This function assumes data has already been converted to long format with an "Abundance" column.
#' @description This function creates an alluvial plot based on input data, which can be
#' either absolute or relative abundance data.
#' @source Built using ggalluvial, ggplot2, and dplyr for visualization of microbial abundance dynamics.
#' @param data A data frame containing abundance and categorical variables.
#' @param axes A character vector specifying the categorical variables for the x-axis.
#' @param abundance_threshold A numeric value specifying the minimum abundance
#' required for an entity to be included in the plot. Default is `10000`.
#' @param fill_variable A string specifying the variable to be used for fill colors. Default is `"Phylum"`.
#' @param silent Logical. If `TRUE`, suppresses warnings. Default is `TRUE`.
#' @param abundance_type A string specifying the type of abundance: `"absolute"` or `"relative"`.
#' @param total_reads Numeric, total number of reads for relative abundance calculation. Default is `NULL`.
#' @param top_taxa Integer. The number of top abundant taxa to retain. Default is `NULL`.
#' @param facet_vars A character vector specifying variables to facet by. Default is `NULL`.
#' @param text_size Numeric, size of text labels. Default is `4`.
#' @param legend_ncol Integer, number of columns for the legend. Default is `1`.
#' @param custom_colors A named vector specifying colors for taxa. Default is `color_palette$MG`.
#' @param color_mapping A named vector of colors for taxa, overriding `custom_colors`. Default is `NULL`.
#' @return A `ggplot2` object representing an alluvial plot.
#' @importFrom ggplot2 ggplot aes geom_text theme scale_x_discrete scale_fill_manual ylab ggtitle guides guide_legend facet_grid
#' @importFrom dplyr group_by summarise mutate ungroup arrange desc pull across filter
#' @importFrom ggalluvial is_alluvia_form geom_alluvium geom_stratum
#' @importFrom magrittr %>%
#' @importFrom grid unit
#' @importFrom rlang sym as_name
#' @examples
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   data("physeq_16SOTU", package = "DspikeIn")
#'
#'   # Convert phyloseq object to long format
#'   pps_Abs <- get_long_format_data(physeq_16SOTU)
#'
#'   # Example of total reads calculation for relative abundance
#'   total_reads <- sum(pps_Abs$Abundance)
#'   print(paste("Total reads:", total_reads))
#'
#'   # Generate an alluvial plot for relative abundance using the extended palette
#'   alluvial_plot_rel <- alluvial_plot(
#'     data = pps_Abs,
#'     axes = c("Env.broad.scale", "Host.genus", "Diet"),
#'     abundance_threshold = 0.01,
#'     fill_variable = "Phylum",
#'     silent = TRUE,
#'     abundance_type = "relative",
#'     top_taxa = 10,
#'     text_size = 4,
#'     legend_ncol = 1,
#'     custom_colors = DspikeIn::color_palette$cool_MG # Extended color palette
#'   )
#'
#'   print(alluvial_plot_rel)
#'
#'   # Convert phyloseq object to TreeSummarizedExperiment (TSE)
#'   tse_data <- convert_phyloseq_to_tse(physeq_16SOTU)
#'   tse_long <- get_long_format_data(tse_data)
#'
#'   # Generate an alluvial plot for absolute abundance
#'   alluvial_plot_abs <- alluvial_plot(
#'     data = tse_long,
#'     axes = c("Env.broad.scale", "Host.genus", "Diet"),
#'     abundance_threshold = 8000,
#'     fill_variable = "Phylum",
#'     silent = TRUE,
#'     abundance_type = "absolute",
#'     top_taxa = 10,
#'     text_size = 4,
#'     legend_ncol = 1,
#'     custom_colors = DspikeIn::color_palette$cool_MG
#'   )
#'
#'   print(alluvial_plot_abs)
#' }
#'
#' @export
alluvial_plot <- function(data, axes = NULL, abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE,
                          abundance_type = "absolute", total_reads = NULL, top_taxa = NULL,
                          facet_vars = NULL, text_size = 4, legend_ncol = 1,
                          custom_colors = color_palette$MG, color_mapping = NULL) {
  # Remove rows with NA values
  data <- stats::na.omit(data)

  # Ensure required columns exist
  if (!is.null(axes) && !all(axes %in% names(data))) {
    stop("Error: Some specified axes are not in the data.")
  }

  if (!fill_variable %in% names(data)) {
    stop("Error: Fill variable is not in the data.")
  }

  # Ensure "Sample" column exists for per-sample normalization
  if (!"Sample" %in% colnames(data)) {
    stop("Error: 'Sample' column is missing. Each row must belong to a specific sample.")
  }

  # Remove rows with NA values in specified axes and abundance
  data <- data[stats::complete.cases(data[, c("Abundance", axes)]), ]

  # Fix for Relative Abundance (Ensuring Proper Normalization)
  if (abundance_type == "relative") {
    if (!is.null(total_reads)) abundance_threshold <- abundance_threshold / total_reads

    # Normalize within each sample first (ensuring each sample sums to 100%)
    data <- data %>%
      dplyr::group_by(Sample) %>%
      dplyr::mutate(RelativeAbundance = Abundance / sum(Abundance) * 100) %>%
      dplyr::ungroup()

    # Summarize within each factor separately
    data <- data %>%
      dplyr::group_by(across(all_of(axes)), !!rlang::sym(fill_variable)) %>%
      dplyr::summarise(RelativeAbundance = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop")

    # Apply filtering
    data <- data[data$RelativeAbundance > abundance_threshold, ]
    abundance_column <- "RelativeAbundance"
  }

  # Fix for Absolute Abundance (Ensure Correct Grouping Across Axes)
  if (abundance_type == "absolute") {
    data <- data %>%
      dplyr::group_by(across(all_of(axes)), !!rlang::sym(fill_variable)) %>%
      dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

    # Apply the threshold filter
    data <- data[data$Abundance > abundance_threshold, ]
    abundance_column <- "Abundance"
  }

  # Fix: Properly Select Top Taxa
  if (!is.null(top_taxa)) {
    top_taxa_names <- data %>%
      dplyr::group_by(!!rlang::sym(fill_variable)) %>%
      dplyr::summarise(TotalAbundance = sum(!!rlang::sym(abundance_column))) %>%
      dplyr::slice_max(n = top_taxa, order_by = TotalAbundance) %>%
      dplyr::pull(!!rlang::sym(fill_variable))

    data <- data %>%
      dplyr::filter(!!rlang::sym(fill_variable) %in% top_taxa_names)
  }

  # Order the levels of the fill variable based on abundance
  data <- data %>%
    dplyr::group_by(!!rlang::sym(fill_variable)) %>%
    dplyr::mutate(TotalAbundance = sum(!!rlang::sym(abundance_column))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(TotalAbundance))
  data[[fill_variable]] <- factor(data[[fill_variable]], levels = unique(data[[fill_variable]]))

  # Remove any NA levels from the fill_variable
  data[[fill_variable]] <- droplevels(data[[fill_variable]])

  # Set colors to use
  if (!is.null(color_mapping)) {
    filtered_color_mapping <- color_mapping[names(color_mapping) %in% unique(data[[fill_variable]])]

    # Warn if any colors are missing for specific taxa
    if (any(is.na(filtered_color_mapping))) {
      warning("\U0000274C There are missing values in the color mapping.")
    }

    color_palette <- filtered_color_mapping

    # Explicitly remove 'NA' from the color mapping
    if ("NA" %in% names(color_palette)) {
      color_palette <- color_palette[!names(color_palette) %in% "NA"]
    }
  } else if (!is.null(custom_colors)) {
    color_palette <- custom_colors
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("\U0000274C Insufficient values in manual scale.")
    }
  } else {
    color_palette <- color_palette$MG # Use MG by default
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("\U0000274C Insufficient values in manual scale.")
    }
  }

  # Create the aes mappings for each axis if axes are specified
  if (!is.null(axes)) {
    axis_mapping <- setNames(lapply(seq_along(axes), function(i) rlang::sym(axes[i])), paste0("axis", seq_along(axes)))

    # Create the plot with axes
    AllE <- ggplot2::ggplot(data, ggplot2::aes(
      y = !!rlang::sym(abundance_column),
      !!!axis_mapping, # Dynamically map the axes
      fill = !!rlang::sym(fill_variable)
    )) +
      ggalluvial::geom_alluvium(width = 0.5, alpha = 0.9, decreasing = TRUE) +
      ggalluvial::geom_stratum(alpha = 0.7, width = 0.4, fill = "gray87", color = "gray50") +
      ggplot2::geom_text(stat = ggalluvial::StatStratum, ggplot2::aes(label = ggplot2::after_stat(stratum)), size = text_size, color = "black") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(color = "black"), # Add X-axis line
        axis.line.y = ggplot2::element_line(color = "black"), # Add Y-axis line
        axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(size = 14, margin = ggplot2::margin(r = 10)),
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      ggplot2::scale_x_discrete(limits = axes, expand = c(.1, .1)) +
      ggplot2::scale_fill_manual(values = color_palette, na.translate = FALSE) + # Ignore NA in legend
      ggplot2::ylab(if (abundance_type == "relative") "Relative Abundance" else "Absolute Abundance") +
      ggplot2::xlab("") + # Remove "Factors" from x-axis when axes are NULL
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
  }

  # Add faceting if specified
  if (!is.null(facet_vars)) {
    AllE <- AllE + ggplot2::facet_grid(stats::reformulate(facet_vars))
  }

  return(AllE)
}

# Example:
# Convert a phyloseq object to a long-format data frame
# pps_Abs <- phyloseq::psmelt(ps)
# tse_long <- get_long_format_data(tse_data) # TSE format
# pps_Abs <- get_long_format_data(SalITSOTU_TSE)

# Example of total reads calculation for relative abundance
# total_reads <- sum(pps_Abs$Abundance)

# Generate an alluvial plot using the extended palette from your package
# alluvial_plot_abs <- alluvial_plot(
#   data = pps_Abs,
#   axes = c("Env.broad.scale", "Host.genus", "Diet"),
#   abundance_threshold = 10000,
#   fill_variable = "Phylum",
#   silent = TRUE,
#   abundance_type = "absolute",
#   top_taxa = 10,
#   text_size = 4,
#   legend_ncol = 1,
#   custom_colors = DspikeIn::color_palette$extended_palette )

# Print the alluvial plot for absolute abundance
# print(alluvial_plot_abs)
