#' Plot Bar Abundance for Phyloseq Data
#'
#' This function generates bar plots for phyloseq data at a specified taxonomic level,
#' with options to customize the appearance/size/legend, relativize the data or plot absolute abundance, facet the plots, and save the plots.
#'
#' @param physeq A phyloseq object containing the microbiome data.
#' @param level A character string specifying the taxonomic level to plot (e.g., "Genus", "Family").
#' @param color A character vector specifying colors to use for the different taxa. Default is NULL, which will use the MG color palette.
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
#' \dontrun{
#' if (interactive()) {
#'   # Plot relativized abundance or select relativize = FALSE to plot absolute abundance
#'   plot <- plotbar_abundance(physeq_16SOTU, 
#'                             level = "Family", 
#'                             group = c("Diet", "Host.species", "Ecoregion.III"), 
#'                             x_axis_var = "Diet", 
#'                             top = 10, x_size = 10, y_size = 10, 
#'                             legend_key_size = 2, legend_text_size = 14,
#'                             legend_nrow = 10, relativize = TRUE, 
#'                             output_prefix = "relativized_abundance_plot", 
#'                             facet_var = "Ecoregion.III", scales = "free_x")
#'   print(plot)
#' }
#' }
#' @importFrom phyloseq psmelt tax_table subset_taxa sample_data otu_table taxa_are_rows
#' @importFrom dplyr group_by summarise arrange pull across
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_y_continuous ylab xlab scale_fill_manual theme_minimal theme element_text element_line element_blank unit guides guide_legend facet_wrap facet_grid aes
#' @importFrom scales percent_format comma
#' @importFrom stats complete.cases as.formula
#' @importFrom utils head
#' @export
plotbar_abundance <- function(physeq, level = "Genus", color = NULL, group = NULL, x_axis_var = NULL, top = 20, return = FALSE, x_size = 12, y_size = 12, legend_key_size = 1.5, legend_text_size = 12, legend_nrow = 20, relativize = TRUE, output_prefix = NULL, facet_var = NULL, scales = "fixed") {
  # Filter out taxa with missing or empty values before melting
  physeq <- phyloseq::subset_taxa(physeq, apply(phyloseq::tax_table(physeq), 1, function(x) all(x != "" & !is.na(x))))
  
  # Melt the phyloseq object into a long-format data frame
  pm <- phyloseq::psmelt(physeq)
  
  # Fix column name conflicts by renaming the taxonomic column if it conflicts with sample data
  tax_levels <- phyloseq::tax_table(physeq)@.Data
  sample_vars <- phyloseq::sample_data(physeq)@names
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
    color <- MG[1:len]  # Use the MG color palette
  }
  
  # Remove rows with NA values in the specified variables
  variables_to_check <- c(group, x_axis_var, facet_var, level)
  pm <- pm[stats::complete.cases(pm[, variables_to_check]), ]
  
  # Summarize abundance by group and taxonomic level
  gv <- dplyr::group_by(pm, dplyr::across(all_of(c(group, level)))) %>%
    dplyr::summarise(summary = sum(Abundance), .groups = 'drop')
  gv <- as.data.frame(gv)
  
  # Select top taxa
  gvy <- dplyr::group_by(pm, dplyr::across(all_of(level))) %>%
    dplyr::summarise(summary = sum(Abundance), .groups = 'drop')
  gvy <- gvy[order(gvy$summary, decreasing = TRUE), ]
  sel <- utils::head(gvy, top) %>% dplyr::pull(level)
  
  # Filter summarized data to include only top taxa
  gv <- gv[gv[[level]] %in% sel, ]
  
  # Create the plot
  p <- ggplot2::ggplot(gv, ggplot2::aes_string(x = x_axis_var[1], y = "summary", fill = level)) +
    ggplot2::geom_bar(stat = "identity", position = if (relativize) "fill" else "stack") +
    ggplot2::scale_y_continuous(labels = if (relativize) scales::percent_format() else scales::comma, expand = c(0, 0.01)) +
    ggplot2::ylab(if (relativize) "Percentage" else "Abundance") +
    ggplot2::xlab("") +
    ggplot2::scale_fill_manual(values = color, name = gsub("tax_", "", level)) +  # Rename legend title
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, size = x_size, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = y_size),
      legend.key.size = ggplot2::unit(legend_key_size, "cm"),  # Adjust unit here to "cm"
      legend.text = ggplot2::element_text(size = legend_text_size, face = "italic"),
      axis.line = ggplot2::element_line(),  # Ensure the axes lines are present
      panel.grid = ggplot2::element_blank(),  # Remove grid lines
      panel.border = ggplot2::element_blank(),  # Remove panel border
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_nrow))
  
  # Add faceting if specified
  if (!is.null(facet_var)) {
    if (length(facet_var) == 1) {
      p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = scales)
    } else if (length(facet_var) == 2) {
      p <- p + ggplot2::facet_grid(stats::as.formula(paste(facet_var[1], "~", facet_var[2])), scales = scales)
    } else {
      stop("facet_var should be a character vector of length 1 or 2.")
    }
  }
  
  # Handle multiple grouping variables for x-axis
  if (length(x_axis_var) > 1) {
    gv$x_axis_combined <- apply(gv[, x_axis_var], 1, paste, collapse = " - ")
    p <- p + ggplot2::aes_string(x = "x_axis_combined")
  } else {
    p <- p + ggplot2::scale_x_discrete(expand = c(0, 0))
  }
  
  # Save plot if output prefix is provided
  if (!is.null(output_prefix)) {
    # Save plot as PNG
    png_filename <- paste0(output_prefix, ".png")
    ggplot2::ggsave(png_filename, plot = p, width = 10, height = 8)
    
    # Save plot as PDF
    pdf_filename <- paste0(output_prefix, ".pdf")
    ggplot2::ggsave(pdf_filename, plot = p, width = 10, height = 8)
  }
  
  # Return the plot or summarized data frame
  if (isTRUE(return)) {
    return(gv)
  } else {
    return(p)
  }
}
# Example usage:
# plot <- plotbar_abundance(
#   physeq = physeq16SOTU,
#   level = "Family",
#   color = color_palette$MG,
#   group = c("Diet", "Host.species", "Host.genus","Ecoregion.III"),
#   x_axis_var = "Host.species",
#   top = 20,
#   x_size = 10,
#   y_size = 10,
#   legend_key_size = 1,
#   legend_text_size = 11,
#   legend_nrow = 20,
#   relativize = T,
#   output_prefix = "rel.abundance_plot",
#   facet_var = "Diet",
#   scales = "free_x"
# )
# print(plot)+my_custom_theme()+ggtitle("across ecoregions")