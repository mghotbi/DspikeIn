#' Plot Bar Abundance for Phyloseq Data
#'
#' This function generates bar plots for phyloseq data at a specified taxonomic level,
#' with options to customize the appearance, relativize the data, and save the plots.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param level A character string specifying the taxonomic level to plot (e.g., "Genus", "Family").
#' @param color A character vector specifying colors to use for the different taxa. Default is NULL.
#' @param group A character string specifying the grouping variable. Default is NULL.
#' @param top An integer specifying the number of top taxa to display. Default is 20.
#' @param return A logical indicating whether to return the summarized data frame. Default is FALSE.
#' @param x_size An integer specifying the size of the x-axis text. Default is 12.
#' @param y_size An integer specifying the size of the y-axis text. Default is 12.
#' @param legend_key_size A numeric value specifying the size of the legend keys. Default is 1.5.
#' @param legend_text_size An integer specifying the size of the legend text. Default is 12.
#' @param legend_nrow An integer specifying the number of rows for the legend. Default is 20.
#' @param relativize A logical indicating whether to relativize the data. Default is TRUE.
#' @param output_prefix A character string specifying the prefix for the output file names. Default is NULL.
#' @return A ggplot object if return is FALSE, otherwise a data frame with summarized data.
#' @examples
#' # Plot relativized abundance or select relativize = FALSE to plot absolute abundance
#' plot <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", 
#' group = "Env_broad_scale.x", top = 10, x_size = 10, y_size = 10, legend_key_size = 2, legend_text_size = 14,
#' legend_nrow = 10, relativize = TRUE, output_prefix = "relativized_abundance_plot")
#' print(plot)
#'
#' # Plot non-relativized (absolute) abundance
#' plot_absolute <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, 
#' level = "Family", group = "Env_broad_scale.x", top = 10, x_size = 10, y_size = 10, 
#' legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE, 
#' output_prefix = "non_relativized_abundance_plot")
#' print(plot_absolute)
#'
#' # Return summarized data
#' data <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, level = "Family", 
#' group = "Env_broad_scale.x", top = 10, return = TRUE)
#' @export
plotbar_abundance <- function(physeq, level = "Genus", color = NULL, group = NULL, top = 20, return = FALSE, x_size = 12, y_size = 12, legend_key_size = 1, legend_text_size = 12, legend_nrow = 20, relativize = TRUE, output_prefix = NULL) {
  # Load necessary libraries
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  
  # Filter out taxa with missing or empty values before melting
  physeq <- subset_taxa(physeq, apply(tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
  
  # Melt the phyloseq object into a long-format data frame
  pm <- psmelt(physeq)
  
  # Generate colors if not provided
  if (is.null(color)) {
    len <- length(unique(pm[, level]))
    color <- rainbow(len)
  }
  
  # Determine grouping variable
  if (is.null(group)) {
    group_var <- c("Sample", level)
  } else {
    group_var <- c(group, level)
  }
  
  # Summarize abundance by group and taxonomic level
  gv <- pm %>% group_by_at(vars(one_of(group_var))) %>% summarise(summary = sum(Abundance))
  gv <- as.data.frame(gv)
  
  # Select top taxa
  gvy <- pm %>% group_by_at(vars(one_of(level))) %>% summarise(summary = sum(Abundance))
  gvy <- gvy[order(gvy$summary, decreasing = TRUE), ]
  sel <- gvy %>% head(top) %>% select(!!sym(level)) %>% pull(1)
  
  # Filter summarized data to include only top taxa
  gv <- gv[gv[, level] %in% sel, ]
  
  # Create the plot
  if (relativize) {
    p <- ggplot(gv, aes_string(x = if (is.null(group)) "Sample" else group, y = "summary", fill = level)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
      ylab("Percentage")
  } else {
    p <- ggplot(gv, aes_string(x = if (is.null(group)) "Sample" else group, y = "summary", fill = level)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(expand = c(0, 0.01)) +
      ylab("Abundance")
  }
  
  # Customize the plot appearance
  p <- p + scale_fill_manual(values = color) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = x_size, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = y_size),
          legend.key.size = unit(legend_key_size, "lines"),
          legend.text = element_text(size = legend_text_size, face = "italic"),
          axis.line = element_line(),  # Ensure the axes lines are present
          panel.grid = element_blank(),  # Remove grid lines
          panel.border = element_blank(),  # Remove panel border
          axis.ticks.x = element_blank()) +
    xlab("") +
    guides(fill = guide_legend(nrow = legend_nrow))
  
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
    return(pm[, c("OTU", "Abundance", group_var)])
  } else {
    return(p)
  }
}

# Example usage:
# Plot relativized abundance
# plot <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, 
#level = "Family", group = "Env_broad_scale.x", top = 10, x_size = 10, y_size = 10,
#legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = TRUE, 
#output_prefix = "relativized_abundance_plot")
# print(plot)

# Plot non-relativized (absolute) abundance
# plot_absolute <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, 
#level = "Family", group = "Env_broad_scale.x", top = 10, x_size = 10, y_size = 10,
#legend_key_size = 2, legend_text_size = 14, legend_nrow = 10, relativize = FALSE,
#output_prefix = "non_relativized_abundance_plot")
# print(plot_absolute)

# Return summarized data
# data <- plotbar_abundance(physeq_16S_adj_scaled_absolute_abundance, 
#level = "Family", group = "Env_broad_scale.x", top = 10, return = TRUE)
