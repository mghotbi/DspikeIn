#' Generate an Alluvial Plot from a Long-Format Data Frame
#'
#' This function generates an alluvial plot using a long-format data frame where microbial abundances are plotted across categorical variables (axes).
#' It supports both relative and absolute abundance types, filters based on abundance thresholds, and allows for displaying the top taxa.
#' Additional customization options include faceting, custom colors, and setting the legend.
#' The function supports displaying multiple axes or plotting without axes when \code{axes = NULL}.
#' 
#' The function uses the \code{ggalluvial} package to construct the plot and relies on \code{ggplot2} for customization.
#' 
#' @param data A data frame containing microbial data, including abundance and categorical variables.
#' @param axes A character vector specifying the axes (columns) to include in the alluvial plot. If \code{NULL}, no axes will be shown on the x-axis.
#' @param abundance_threshold A numeric value specifying the minimum abundance threshold for including samples. Default is 10000 for absolute abundance and a relative threshold for relative abundance (e.g., 0.001 for 0.1\%).
#' @param fill_variable A character string specifying the variable to use for the fill color in the alluvial plot. Default is "Phylum".
#' @param silent A logical indicating whether to suppress messages from \code{ggalluvial::is_alluvia_form}. Default is \code{TRUE}.
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "absolute".
#' @param total_reads A numeric value specifying the total number of reads for relative abundance calculation. Default is \code{NULL}.
#' @param top_taxa An integer specifying the number of top taxa to display. Default is \code{NULL}, meaning all taxa will be displayed.
#' @param facet_vars A character vector specifying variables for faceting. Default is \code{NULL} (no faceting).
#' @param text_size A numeric value specifying the size of the text labels inside the plot. Default is 4.
#' @param legend_ncol An integer specifying the number of columns for the legend. Default is 1.
#' @param custom_colors A character vector specifying custom colors for the fill variable. Can use \code{color_palette$MG} or \code{color_palette$extended_palette}. Default is \code{color_palette$MG}.
#' @param color_mapping A named character vector specifying specific colors for taxa. Default is \code{NULL}.
#' @return A \code{ggplot2} object representing the alluvial plot.
#' @importFrom ggplot2 ggplot aes geom_text theme scale_x_discrete scale_fill_manual ylab ggtitle guides guide_legend facet_grid
#' @importFrom dplyr group_by summarise mutate ungroup arrange desc pull across filter
#' @importFrom ggalluvial is_alluvia_form geom_alluvium geom_stratum
#' @importFrom magrittr %>%
#' @importFrom grid unit
#' @examples
#' \dontrun{
#' # Example data
#' data(physeq_16SOTU)  # Replace with actual data loading
#' pps_Abs <- phyloseq::psmelt(physeq_16SOTU)
#'
#' # Example of total reads calculation for relative abundance
#' total_reads <- sum(pps_Abs$Abundance)
#'
#' # Generate an alluvial plot using the extended palette
#' alluvial_plot_abs <- alluvial_plot(
#'   data = pps_Abs, 
#'   axes = c("Env.broad.scale", "Host.genus", "Diet"), 
#'   abundance_threshold = 10000, 
#'   fill_variable = "Phylum", 
#'   silent = TRUE, 
#'   abundance_type = "absolute", 
#'   top_taxa = 10, 
#'   text_size = 4, 
#'   legend_ncol = 1, 
#'   custom_colors = color_palette$extended_palette  # Use the extended palette from your package
#' )
#'
#' # Print the alluvial plot for absolute abundance
#' print(alluvial_plot_abs)
#' }
#' @export
alluvial_plot <- function(data, axes = NULL, abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE, 
                          abundance_type = "absolute", total_reads = NULL, top_taxa = NULL, 
                          facet_vars = NULL, text_size = 4, legend_ncol = 1, 
                          custom_colors = color_palette$MG, color_mapping = NULL) {
  
  # Remove rows with NA values from the entire dataset
  data <- stats::na.omit(data)
  
  # Ensure axes exist in the data if provided
  if (!is.null(axes) && !all(axes %in% names(data))) {
    stop("Some specified axes are not present in the data.")
  }
  
  # Ensure the fill variable exists in the data
  if (!fill_variable %in% names(data)) {
    stop("Fill variable is not present in the data.")
  }
  
  # Convert to relative threshold if specified
  if (abundance_type == "relative" && !is.null(total_reads)) {
    abundance_threshold <- abundance_threshold / total_reads
  }
  
  # Remove rows with NA values in the specified axes and abundance
  data <- data[stats::complete.cases(data[, c("Abundance", axes)]), ]
  
  # Filter out samples with abundance below the threshold
  if (abundance_type == "relative") {
    data <- data %>%
      dplyr::group_by(dplyr::across(all_of(axes))) %>%
      dplyr::mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
      dplyr::ungroup()
    data <- data[data$RelativeAbundance > abundance_threshold, ]
    abundance_column <- "RelativeAbundance"
  } else {
    data <- data[data$Abundance > abundance_threshold, ]
    abundance_column <- "Abundance"
  }
  
  # Filter top taxa if specified
  if (!is.null(top_taxa)) {
    top_taxa_names <- data %>% 
      dplyr::group_by(!!rlang::sym(fill_variable)) %>%
      dplyr::summarise(TotalAbundance = sum(!!rlang::sym(abundance_column))) %>%
      dplyr::top_n(n = top_taxa, wt = TotalAbundance) %>%
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
  
  # Validate the data structure using ggalluvial::is_alluvia_form if not silent
  if (!silent) {
    ggalluvial::is_alluvia_form(as.data.frame(data), axes = axes, silent = TRUE)
  }
  
  # Set colors to use
  if (!is.null(color_mapping)) {
    filtered_color_mapping <- color_mapping[names(color_mapping) %in% unique(data[[fill_variable]])]
    
    # Warn if any colors are missing for specific taxa
    if (any(is.na(filtered_color_mapping))) {
      warning("There are missing values in the color mapping.")
    }
    
    color_palette <- filtered_color_mapping
    
    # Explicitly remove 'NA' from the color mapping
    if ("NA" %in% names(color_palette)) {
      color_palette <- color_palette[!names(color_palette) %in% "NA"]
    }
    
  } else if (!is.null(custom_colors)) {
    color_palette <- custom_colors
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("Insufficient values in manual scale.")
    }
  } else {
    color_palette <- color_palette$MG  # Use MG by default
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("Insufficient values in manual scale.")
    }
  }
  
  # Create the aes mappings for each axis if axes are specified
  if (!is.null(axes)) {
    axis_mapping <- setNames(lapply(seq_along(axes), function(i) rlang::sym(axes[i])), paste0("axis", seq_along(axes)))
    
    # Create the plot with axes
    AllE <- ggplot2::ggplot(data, ggplot2::aes(
      y = !!rlang::sym(abundance_column),
      !!!axis_mapping,  # Dynamically map the axes
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
        axis.line.x = ggplot2::element_line(color = "black"),  # Add X-axis line
        axis.line.y = ggplot2::element_line(color = "black"),  # Add Y-axis line
        axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(size = 14, margin = ggplot2::margin(r = 10)),
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      ggplot2::scale_x_discrete(limits = axes, expand = c(.1, .1)) +
      ggplot2::scale_fill_manual(values = color_palette, na.translate = FALSE) +  # Ignore NA in legend
      ggplot2::ylab(if (abundance_type == "relative") "Relative Abundance (%)" else "Absolute Abundance") +
      ggplot2::xlab("") +  # Remove "Factors" from x-axis when axes are NULL
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
    
  } else {
    # Create the plot without specific axes
    AllE <- ggplot2::ggplot(data, ggplot2::aes(
      y = !!rlang::sym(abundance_column),
      fill = !!rlang::sym(fill_variable)
    )) +
      ggalluvial::geom_alluvium(width = 0.5, alpha = 0.9, decreasing = TRUE) +
      ggalluvial::geom_stratum(alpha = 0.7, width = 0.4, fill = "gray80", color = "gray50") +
      ggplot2::geom_text(stat = ggalluvial::StatStratum, ggplot2::aes(label = ggplot2::after_stat(stratum)), size = text_size, color = "black") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right", 
        legend.title = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 10),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(color = "black"),  # Add X-axis line
        axis.line.y = ggplot2::element_line(color = "black"),  # Add Y-axis line
        axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(size = 14, margin = ggplot2::margin(r = 10)),
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)
      ) +
      ggplot2::scale_fill_manual(values = color_palette, na.translate = FALSE) +  # Ignore NA in legend
      ggplot2::ylab(if (abundance_type == "relative") "Relative Abundance (%)" else "Abundance") +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
  }
  
  # Add faceting if specified
  if (!is.null(facet_vars)) {
    AllE <- AllE + ggplot2::facet_grid(stats::reformulate(facet_vars))
  }
  
  return(AllE)
}


# Example:
# Load necessary libraries
# library(phyloseq)
# library(ggplot2)
# library(dplyr)
# library(ggalluvial)

# Assuming `color_palette` is already defined in your package with MG and extended_palette

# Convert a phyloseq object to a long-format data frame
# pps_Abs <- phyloseq::psmelt(physeq_16SOTU)

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
#   custom_colors = color_palette$extended_palette  # Use the extended palette from your package
# )

# Print the alluvial plot for absolute abundance
# print(alluvial_plot_abs)
