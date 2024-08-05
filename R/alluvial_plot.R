#' Generate an Alluvial Plot after converting phyloseq object to the long-format dataframe
#'
#' This function creates an alluvial plot to visualize the abundance of common taxa across multiple categorical axes/which reveals factors impact on abundance of taxa.
#' It filters out samples with abundance below the specified threshold and generates the plot using `ggplot2`.
#'
#' @param data A data frame containing the microbial data, including abundance and categorical variables.
#' @param axes A character vector specifying the axes (columns) to include in the alluvial plot.
#' @param abundance_threshold A numeric value specifying the minimum abundance threshold for including samples. Default is 10000.
#' @param fill_variable A character string specifying the variable to use for fill in the alluvial plot. Default is "Phylum".
#' @param silent A logical indicating whether to suppress messages from `is_alluvia_form`. Default is TRUE.
#' @param abundance_type A character string specifying whether to plot "relative" or "absolute" abundance. Default is "absolute".
#' @param total_reads A numeric value specifying the total number of reads for the sample when using relative abundance. Default is NULL.
#' @param top_taxa An integer specifying the number of top taxa to display. Default is NULL, which means all taxa will be displayed.
#' @param facet_vars A character vector specifying the variables for faceting. Default is NULL.
#' @param text_size A numeric value specifying the size of the text inside the plot. Default is 4.
#' @param legend_ncol An integer specifying the number of columns for the legend. Default is 1.
#' @return A ggplot2 object representing the alluvial plot.
#' @examples
#' # Load necessary libraries
#' library(phyloseq)
#' library(ggplot2)
#' library(dplyr)
#' library(ggalluvial)
#' 
#' # Assuming 'ps' is your phyloseq object
#' ps <- phyloseq object
#' 
#' # Melt the phyloseq object into a long-format data frame
#' pps <- psmelt(ps)
#' 
#' # Define total reads for relative abundance calculation
#' total_reads <- sum(pps$Abundance)  # Calculate total reads from the data
#' 
#' # Generate alluvial plot for absolute abundance
#' alluvial_plot_abs <- alluvial_plot(data = pps, axes = c("Animal.ecomode", "Host.species", "Result"), 
#'                                    abundance_threshold = 10000, fill_variable = "Phylum", 
#'                                    silent = TRUE, abundance_type = "absolute", top_taxa = 10,
#'                                    text_size = 4, legend_ncol = 1)
#' print(alluvial_plot_abs)
#' 
#' # Generate alluvial plot for relative abundance
#' alluvial_plot_rel <- alluvial_plot(data = pps, axes = c("Animal.ecomode", "Host.species", "Result"), 
#'                                    abundance_threshold = 10000, fill_variable = "Phylum", 
#'                                    silent = TRUE, abundance_type = "relative", 
#'                                    total_reads = total_reads, top_taxa = 10, text_size = 4, legend_ncol = 1)
#' print(alluvial_plot_rel)
#' @export
alluvial_plot <- function(data, axes, abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE, abundance_type = "absolute", total_reads = NULL, top_taxa = NULL, facet_vars = NULL, text_size = 4, legend_ncol = 1) {
  # Load necessary libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required but not installed.")
  }
  
  library(ggplot2)
  library(dplyr)
  library(ggalluvial)
  
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
  
  # Check if is_alluvia_form needs to be called
  if (!silent) {
    is_alluvia_form(as.data.frame(data), axes = axes, silent = silent)
  }
  
  # Get unique fill values
  unique_fill_values <- unique(data[[fill_variable]])
  num_unique_fill_values <- length(unique_fill_values)
  
  # Check if MG has enough colors
  if (num_unique_fill_values > length(MG)) {
    stop(paste("Insufficient values in manual scale. ", num_unique_fill_values, " needed but only ", length(MG), " provided.", sep = ""))
  }
  
  # Create the alluvial plot
  AllE <- ggplot(data, aes(y = .data[[abundance_column]], !!!setNames(lapply(axes, as.name), paste0("axis", seq_along(axes))))) +
    geom_alluvium(aes_string(fill = fill_variable), width = 0.5, alpha = 0.8, decreasing = TRUE) +
    geom_stratum(alpha = 0.5, width = 0.3, fill = "gray80", color = "gray30") +
    geom_label(stat = "stratum", size = text_size, hjust = 0.5, vjust = 0.5, aes_string(label = "after_stat(stratum)"), reverse = FALSE) +
    theme(legend.position = "right") +
    scale_x_discrete(limits = axes, expand = c(.0, .0)) +
    scale_fill_manual(values = MG[1:num_unique_fill_values]) +
    ylab(if (abundance_type == "relative") "Relative Abundance" else "Abundance") +
    ggtitle("Abundance across factors") +
    my_custom_theme() +
    guides(fill = guide_legend(ncol = legend_ncol))
  
  # Add faceting if specified
  if (!is.null(facet_vars)) {
    facet_formula <- as.formula(paste("~", paste(facet_vars, collapse = "+")))
    AllE <- AllE + facet_grid(facet_formula)
  }
  
  return(AllE)
}

#' Custom GGplot2 Theme and Color/Shape Vectors
#'
#' This function generates a custom ggplot2 theme with specific styling options.
#' Additionally, it provides predefined color and shape vectors for consistent plot styling.
#'
#' @return A ggplot2 theme object.
#' @examples
#' # Apply the custom theme to a ggplot
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + 
#'   geom_point(size = 3) + 
#'   my_custom_theme()
#' print(p)
#' @export
my_custom_theme <- function() {
  theme_bw() + 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = 'white', color = "#e1deda"),
      panel.border = element_blank(),
      axis.line.x = element_line(colour = 'black', size = 0.6),
      axis.line.y = element_line(colour = 'black', size = 0.6),
      axis.ticks = element_line(colour = 'black', size = 0.35),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 11, color = "black", face = "bold"), 
      legend.key.size = unit(0.5, 'cm'),
      axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), 
      axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), 
      axis.text.x = element_text(family = "Times New Roman", size = 12, angle = 0, color = "black", face = "bold"), 
      axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      plot.title = element_text(color = "black", size = 12, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
}

#' Predefined color and shape vectors for consistent plot styling
#' @export
MG <- c("#FFFF33", "#FF7F00", "#E41A1C", "firebrick4", "#2e4057", "#984EA3", "#377EB8", "olivedrab3", "#4DAF4A", "#336633", "grey80", "#BB650B", "gold", "#559999", "#7570b3", "#E78AC3", "#A6D854", "#66a61e", "#e6ab02", "#a6761d", "#663300", "#66C2A5", "#0e669b", "#00798c", "dodgerblue4", "steelblue2", "#00AFBB", "#E7B800", "#FC4E07", "lightskyblue4", "green", "red", "#FFF000", "#0099CC", "#FF9933", "#CC9900", "chartreuse1", "#FF3399", "#00FFFF", "#0000CC", "#A37F6F", "#9183E6", "#00AD9A", "#990033", "#909800", "#00FF00", "#17b5b4", "#AED1D6", "#b1010c", "firebrick2", "blue", "navy", "yellow", "brown", "black", "purple", "darkred", "darkgreen", "#82cfd0", "#b2e0e4", "honeydew3", "#8d96a3", "lavender", "#CC6686", "lavenderblush2", "mistyrose3", "#e1deda", "darkgoldenrod", "burlywood", "papayawhip", "wheat4", "cornsilk3", "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream", "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue")

#' @export
MG_shape <- c(19, 3, 1, 2, 9, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 23, 20, 22)

# Example usage
# use phyloseq object
# pss_abs<-psmelt(ps.abs)
# pps_rel<-psmelt(ps.rel)
# # Define total reads for relative abundance calculation
# total_reads <- sum(pps$Abundance)  # Calculate total reads from the data
# # Generate alluvial plot for absolute abundance
# alluvial_plot_rel <- alluvial_plot(data = pps_rel,
#                                    axes = c("Clade.Order","Ecoregion.III","Diet", "Habitat","Animal.ecomode" ),
#                                    abundance_threshold = 3000, fill_variable = "Family",
#                                    silent = TRUE, abundance_type = "absolute",
#                                    top_taxa = 30, text_size = 3, legend_ncol = 1)
# 
# # Generate alluvial plot for relative abundance
# alluvial_plot_abs <- alluvial_plot(data = pss_abs,
#                                    axes = c("Diet", "Host.genus"),
#                                    abundance_threshold = 1000, fill_variable = "Family",
#                                    silent = TRUE, abundance_type = "relative",
#                                    total_reads = total_reads, top_taxa = 10,legend_ncol = 1,
#                                    facet_vars = c("Ecoregion.III"))
# print(alluvial_plot_rel)
