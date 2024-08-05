#' Plot Bar Abundance for Phyloseq Data
#'
#' This function generates bar plots for phyloseq data at a specified taxonomic level,
#' with options to customize the appearance/size/legend, relativize the data or plot absolute abundance, facet the plots, and save the plots.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param level A character string specifying the taxonomic level to plot (e.g., "Genus", "Family").
#' @param color A character vector specifying colors to use for the different taxa. Default is NULL.
#' @param group A character vector specifying the grouping variables. Default is NULL.
#' @param x_axis_var A character vector specifying the variable(s) to be shown on the x-axis. Default is NULL.
#' @param top An integer specifying the number of top taxa to display. Default is 20.
#' @param return A logical indicating whether to return the summarized data frame. Default is FALSE.
#' @param x_size An integer specifying the size of the x-axis text. Default is 12.
#' @param y_size An integer specifying the size of the y-axis text. Default is 12.
#' @param legend_key_size A numeric value specifying the size of the legend keys. Default is 1.5.
#' @param legend_text_size An integer specifying the size of the legend text. Default is 12.
#' @param legend_nrow An integer specifying the number of rows for the legend. Default is 20.
#' @param relativize A logical indicating whether to relativize the data. Default is TRUE.
#' @param output_prefix A character string specifying the prefix for the output file names. Default is NULL.
#' @param facet_var A character vector specifying the faceting variables. Default is NULL.
#' @param scales A character string specifying scales in facet. Options are "fixed", "free", "free_x", "free_y". Default is "fixed".
#' @return A ggplot object if return is FALSE, otherwise a data frame with summarized data.
#' @examples
#' # Plot relativized abundance or select relativize = FALSE to plot absolute abundance
#' plot <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", 
#' group = c("Diet", "Host.species","Ecoregion.III"), x_axis_var = "Diet", 
#' top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14,
#' legend_nrow = 10, relativize = TRUE, output_prefix = "relativized_abundance_plot", 
#' facet_var = "Ecoregion.III", scales = "free_x")
#' print(plot)
#'
#' # Plot non-relativized (absolute) abundance
#' plot_absolute <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, 
#' level = "Family", group = c("Diet", "Host.species","Ecoregion.III"), 
#' x_axis_var = "Env_broad_scale.x", top = 10, x_size = 10, y_size = 10, 
#' legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE, 
#' output_prefix = "non_relativized_abundance_plot", facet_var = "Ecoregion.III", scales = "free_x")
#' print(plot_absolute)
#'
#' # Return summarized data
#' data <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", 
#' group = c("Diet", "Host.species","Ecoregion.III"), 
#' x_axis_var = "Ecoregion.III", top = 10, return = TRUE)
#' @export
plotbar_abundance <- function(physeq, level = "Genus", color = NULL, group = NULL, x_axis_var = NULL, top = 20, return = FALSE, x_size = 12, y_size = 12, legend_key_size = 1.5, legend_text_size = 12, legend_nrow = 20, relativize = TRUE, output_prefix = NULL, facet_var = NULL, scales = "fixed") {
  # Load necessary libraries
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  
  # Filter out taxa with missing or empty values before melting
  physeq <- subset_taxa(physeq, apply(tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
  
  # Melt the phyloseq object into a long-format data frame
  pm <- psmelt(physeq)
  
  # Fix column name conflicts by renaming the taxonomic column if it conflicts with sample data
  tax_levels <- tax_table(physeq)@.Data
  sample_vars <- sample_data(physeq)@names
  conflicting_names <- intersect(colnames(tax_levels), sample_vars)
  if (length(conflicting_names) > 0) {
    for (name in conflicting_names) {
      colnames(pm)[colnames(pm) == name] <- paste0("tax_", name)
    }
    level <- paste0("tax_", level)
  }
  
  # Generate colors if not provided
  if (is.null(color)) {
    len <- length(unique(pm[[level]]))
    color <- MG()[1:len]  # Use the MG color palette function
  }
  
  # Remove rows with NA values in the specified variables
  variables_to_check <- c(group, x_axis_var, facet_var, level)
  pm <- pm[complete.cases(pm[, variables_to_check]), ]
  
  # Summarize abundance by group and taxonomic level
  gv <- pm %>%
    group_by(across(all_of(c(group, level)))) %>%
    summarise(summary = sum(Abundance), .groups = 'drop')
  gv <- as.data.frame(gv)
  
  # Select top taxa
  gvy <- pm %>%
    group_by(across(all_of(level))) %>%
    summarise(summary = sum(Abundance), .groups = 'drop')
  gvy <- gvy[order(gvy$summary, decreasing = TRUE), ]
  sel <- gvy %>% head(top) %>% pull(level)
  
  # Filter summarized data to include only top taxa
  gv <- gv[gv[[level]] %in% sel, ]
  
  # Create the plot
  p <- ggplot(gv, aes(x = !!sym(x_axis_var[1]), y = summary, fill = !!sym(level))) +
    geom_bar(stat = "identity", position = if (relativize) "fill" else "stack") +
    scale_y_continuous(labels = if (relativize) scales::percent_format() else scales::comma, expand = c(0, 0.01)) +
    ylab(if (relativize) "Percentage" else "Abundance") +
    xlab("") +
    scale_fill_manual(values = color, name = gsub("tax_", "", level)) +  # Rename legend title
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, size = x_size, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = y_size),
      legend.key.size = unit(legend_key_size, "cm"),  # Adjust unit here to "cm"
      legend.text = element_text(size = legend_text_size, face = "italic"),
      axis.line = element_line(),  # Ensure the axes lines are present
      panel.grid = element_blank(),  # Remove grid lines
      panel.border = element_blank(),  # Remove panel border
      axis.ticks.x = element_blank()
    ) +
    guides(fill = guide_legend(nrow = legend_nrow))
  
  # Add faceting if specified
  if (!is.null(facet_var)) {
    if (length(facet_var) == 1) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = scales)
    } else if (length(facet_var) == 2) {
      p <- p + facet_grid(as.formula(paste(facet_var[1], "~", facet_var[2])), scales = scales)
    } else {
      stop("facet_var should be a character vector of length 1 or 2.")
    }
  }
  
  # Handle multiple grouping variables for x-axis
  if (length(x_axis_var) > 1) {
    gv$x_axis_combined <- apply(gv[, x_axis_var], 1, paste, collapse = " - ")
    p <- p + aes_string(x = "x_axis_combined")
  } else {
    p <- p + scale_x_discrete(expand = c(0, 0))
  }
  
  # Save plot if output prefix is provided
  if (!is.null(output_prefix)) {
    # Save plot as PNG
    png_filename <- paste0(output_prefix, ".png")
    ggsave(png_filename, plot = p, width = 10, height = 8)
    
    # Save plot as PDF
    pdf_filename <- paste0(output_prefix, ".pdf")
    ggsave(pdf_filename, plot = p, width = 10, height = 8)
  }
  
  # Return the plot or summarized data frame
  if (isTRUE(return)) {
    return(gv)
  } else {
    return(p)
  }
}

# Example usage:

# Salamander_absolute_NospikeSp_pl <- subset_samples(Salamander_absolute_NospikeSp, 
# Diet == "Insectivore" & Host.genus == "Plethodon")
# Salamander_rel_NospikeSp_pl <- subset_samples(Salamander_relative_NospikeSp, 
# Diet == "Insectivore" & Host.genus == "Plethodon")

# Plot non-relativized (absolute) abundance
# plot_absolute <- plotbar_abundance(Salamander_absolute_NospikeSp_pl,
#                                    level = "Family", group = c("Host.species", "Diet","Ecoregion.III","Animal.ecomode"),
#                                    x_axis_var = "Host.species", top = 25, x_size = 10, y_size = 10,
#                                    legend_key_size = 0.5, legend_text_size = 10, legend_nrow = 25, relativize = FALSE,
#                                    output_prefix = "non_relativized_abundance_plot",
#                                    facet_var = c("Ecoregion.III"), scales = "free_x")
# print(plot_absolute)
# 
# # Plot relative abundance
# plot_rel <- plotbar_abundance(Salamander_rel_NospikeSp_pl,
#                               level = "Family", group = c("Host.species", "Diet","Ecoregion.III","Animal.ecomode"),
#                               x_axis_var = "Host.species", top = 25, x_size = 10, y_size = 10,
#                               legend_key_size = 0.5, legend_text_size = 12, legend_nrow = 25, relativize = TRUE,
#                               output_prefix = "relativized_abundance_plot",
#                               facet_var = c("Ecoregion.III"), scales = "free_x")
# print(plot_rel)
