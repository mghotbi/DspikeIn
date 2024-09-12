#' Generate an Alluvial Plot from a Long-Format Dataframe
#'
#' This function generates an alluvial plot using a long-format data frame where microbial abundances are plotted across categorical variables (axes).
#' It supports both relative and absolute abundance types, filters based on abundance thresholds, and allows for displaying the top taxa.
#' Additional customization options include faceting, custom colors, and setting the legend.
#'
#' @param data A data frame containing microbial data, including abundance and categorical variables.
#' @param axes A character vector specifying the axes (columns) to include in the alluvial plot.
#' @param abundance_threshold A numeric value specifying the minimum abundance threshold for including samples. Default is 10000 for absolute abundance and a relative threshold for relative abundance (e.g., 0.001 for 0.1\%).
#' @param fill_variable A character string specifying the variable to use for the fill color in the alluvial plot. Default is "Phylum".
#' @param silent A logical indicating whether to suppress messages from `is_alluvia_form`. Default is TRUE.
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "absolute".
#' @param total_reads A numeric value specifying the total number of reads for relative abundance calculation. Default is NULL.
#' @param top_taxa An integer specifying the number of top taxa to display. Default is NULL, meaning all taxa will be displayed.
#' @param facet_vars A character vector specifying variables for faceting. Default is NULL (no faceting).
#' @param text_size A numeric value specifying the size of the text labels inside the plot. Default is 4.
#' @param legend_ncol An integer specifying the number of columns for the legend. Default is 1.
#' @param custom_colors A character vector specifying custom colors for the fill variable. Default is NULL (use predefined MG colors).
#' @param color_mapping A named character vector specifying specific colors for taxa. Default is NULL.
#'
#' @return A ggplot2 object representing the alluvial plot.
#'
#' @examples
#' library(phyloseq)
#' library(ggplot2)
#' library(dplyr)
#' library(ggalluvial)
#'
#' # Convert a phyloseq object to a long-format data frame
#' pps <- psmelt(physeq_16SOTU)
#'
#' # Calculate total reads for relative abundance
#' total_reads <- sum(pps$Abundance)
#'
#' # Generate alluvial plot for absolute abundance
#' alluvial_plot_abs <- alluvial_plot(
#'   data = pps, 
#'   axes = c("Animal.ecomode", "Host.species", "Result"), 
#'   abundance_threshold = 10000, 
#'   fill_variable = "Phylum", 
#'   silent = TRUE, 
#'   abundance_type = "absolute", 
#'   top_taxa = 10,
#'   text_size = 4, 
#'   legend_ncol = 1
#' )
#' print(alluvial_plot_abs)
#'
#' # Generate alluvial plot for relative abundance
#' alluvial_plot_rel <- alluvial_plot(
#'   data = pps, 
#'   axes = c("Animal.ecomode", "Host.species", "Diet"), 
#'   abundance_threshold = 0.001, 
#'   fill_variable = "Phylum", 
#'   silent = TRUE, 
#'   abundance_type = "relative", 
#'   total_reads = total_reads, 
#'   top_taxa = 10, 
#'   text_size = 4, 
#'   legend_ncol = 1
#' )
#' print(alluvial_plot_rel)
#'
#' @importFrom ggplot2 ggplot aes geom_label theme scale_x_discrete scale_fill_manual ylab ggtitle guides guide_legend facet_grid
#' @importFrom dplyr group_by summarise mutate ungroup arrange desc pull across filter
#' @importFrom ggalluvial is_alluvia_form geom_alluvium geom_stratum
#' @importFrom magrittr %>%
#' @importFrom grid unit
#' @export
alluvial_plot <- function(data, axes, abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE, 
                          abundance_type = "absolute", total_reads = NULL, top_taxa = NULL, 
                          facet_vars = NULL, text_size = 4, legend_ncol = 1, 
                          custom_colors = NULL, color_mapping = NULL) {
  # Ensure specified axes exist in the data
  if (!all(axes %in% names(data))) {
    stop("Some specified axes are not present in the data.")
  }
  
  # Convert to relative threshold if specified
  if (abundance_type == "relative" && !is.null(total_reads)) {
    abundance_threshold <- abundance_threshold / total_reads
  }
  
  # Remove rows with NA values in the specified axes and abundance
  data <- data[complete.cases(data[, c("Abundance", axes)]), ]
  
  # Filter out samples with abundance below the threshold
  if (abundance_type == "relative") {
    data <- data %>%
      group_by(across(all_of(axes))) %>%
      mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
      ungroup()
    data <- data[data$RelativeAbundance > abundance_threshold, ]
    abundance_column <- "RelativeAbundance"
  } else {
    data <- data[data$Abundance > abundance_threshold, ]
    abundance_column <- "Abundance"
  }
  
  # Filter top taxa if specified
  if (!is.null(top_taxa)) {
    top_taxa_names <- data %>%
      group_by(!!sym(fill_variable)) %>%
      summarise(TotalAbundance = sum(!!sym(abundance_column))) %>%
      top_n(n = top_taxa, wt = TotalAbundance) %>%
      pull(!!sym(fill_variable))
    
    data <- data %>%
      filter(!!sym(fill_variable) %in% top_taxa_names)
  }
  
  # Order the levels of the fill variable based on abundance
  data <- data %>%
    group_by(!!sym(fill_variable)) %>%
    mutate(TotalAbundance = sum(!!sym(abundance_column))) %>%
    ungroup() %>%
    arrange(desc(TotalAbundance))
  data[[fill_variable]] <- factor(data[[fill_variable]], levels = unique(data[[fill_variable]]))
  
  # Check if is_alluvia_form needs to be called
  if (!silent) {
    is_alluvia_form(as.data.frame(data), axes = axes, silent = silent)
  }
  
  # Set colors to use
  if (!is.null(color_mapping)) {
    filtered_color_mapping <- color_mapping[names(color_mapping) %in% unique(data[[fill_variable]])]
    color_palette <- filtered_color_mapping
  } else if (!is.null(custom_colors)) {
    color_palette <- custom_colors
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("Insufficient values in manual scale.")
    }
  } else {
    color_palette <- MG
    if (length(unique(data[[fill_variable]])) > length(color_palette)) {
      stop("Insufficient values in manual scale.")
    }
  }
  
  # Create the alluvial plot
  AllE <- ggplot(data, aes(y = .data[[abundance_column]], !!!setNames(lapply(axes, as.name), paste0("axis", seq_along(axes))))) +
    geom_alluvium(aes(fill = .data[[fill_variable]]), width = 0.5, alpha = 0.8, decreasing = TRUE) +
    geom_stratum(alpha = 0.5, width = 0.3, fill = "gray80", color = "gray30") +
    geom_label(stat = "stratum", size = text_size, aes(label = after_stat(stratum)), reverse = FALSE) +
    theme(legend.position = "right") +
    scale_x_discrete(limits = axes, expand = c(.0, .0)) +
    scale_fill_manual(values = color_palette) +
    ylab(if (abundance_type == "relative") "Relative Abundance (%)" else "Abundance") +
    ggtitle("Abundance across factors") +
    my_custom_theme() +
    guides(fill = guide_legend(ncol = legend_ncol))
  
  # Add faceting if specified
  if (!is.null(facet_vars)) {
    AllE <- AllE + facet_grid(reformulate(facet_vars))
  }
  
  return(AllE)
}

# Example usage
# pps_rel<-psmelt(physeq_16SOTU)
# total_reads <- sum(pps_rel$Abundance)  # Calculate total reads from the data
# MG<-color_palette$MG
# # Generate alluvial plot for relative abundance
# alluvial_plot_rel <- alluvial_plot(
# data = pps_rel, 
# axes = c("Animal.ecomode", "Host.species", "Diet"), 
# abundance_threshold = 0.001, 
# fill_variable = "Phylum", 
# silent = TRUE, 
# abundance_type = "relative", 
# total_reads = total_reads, 
# top_taxa = 10, 
# text_size = 4, 
# legend_ncol = 1)
# print(alluvial_plot_rel)
