#' Create an Alluvial Plot
#'
#' This function generates an alluvial plot from the given dataset, illustrating the abundance
#' across specified factors. Samples with abundance below the given threshold are filtered out.
#'
#' @param data A dataframe containing the data to be plotted.
#' @param axes A character vector specifying the columns of the data to be used as axes in the alluvial plot. Default is c("Host_genus", "Ecoregion_III", "Diet").
#' @param abundance_threshold A numeric value specifying the minimum abundance threshold for samples to be included in the plot. Default is 10000.
#' @param fill_variable A character string specifying the variable used for filling the alluvia. Default is "Phylum".
#' @param silent A logical value indicating whether to suppress warnings and messages. Default is TRUE.
#' @return A ggplot object representing the alluvial plot.
#' @examples
#' # `meli` is a preprocessed dataframe as described in the example
#' alluvial_plot <- alluvial_plot(data = meli, axes = c("Abundance", "Diet", "Host_genus", "Ecoregion_III"), abundance_threshold = 20000, silent = TRUE)
#' # Display the plot
#' print(alluvial_plot)
#' @export
alluvial_plot <- function(data, axes = c("Host_genus", "Ecoregion_III", "Diet"), abundance_threshold = 10000, fill_variable = "Phylum", silent = TRUE) {
  # Check if input data is a dataframe
  if (!is.data.frame(data)) {
    stop("Error: data must be a dataframe.")
  }
  
  # Check if the specified axes are in the dataframe
  if (!all(axes %in% colnames(data))) {
    stop("Error: Not all specified axes are present in the dataframe.")
  }
  
  # Filter out samples with abundance below the threshold
  data <- data[data$Abundance > abundance_threshold, ]
  
  # Check if is_alluvia_form needs to be called
  if (!silent) {
    is_alluvia_form(as.data.frame(data), axes = axes, silent = silent)
  }
  
  # Create the alluvial plot dynamically using the specified axes
  AllE <- ggplot(data, aes_string(y = "Abundance", axis1 = axes[1], axis2 = axes[2], axis3 = axes[3])) +
    geom_alluvium(aes_string(fill = fill_variable), width = 1/5, alpha = 0.9, decreasing = TRUE) +
    geom_stratum(alpha = 0.5, width = 1/5, fill = "gray", color = "black") +
    geom_label(stat = "stratum", size = 4, aes(label = after_stat(stratum)), reverse = FALSE) +
    theme(legend.position = "none") +
    scale_x_discrete(limits = axes, expand = c(0, 0)) +
    scale_fill_manual(values = MG) +  # Use extended_palette instead of MG if you have many items
    ylab("Abundance") +
    ggtitle("Abundance across factors") +
    mytheme
  
  return(AllE)
}


# Example usage:
#SF <- phyloseq::subset_samples(physeqASV16, ecoregion.III == "Interior Plateau" & env.broad.scale=="wild")
#AFN <- subset_samples(SF, !is.na(ecoregion.III) & !is.na(env.broad.scale)  & !is.na(host.genus) & !is.na(animal.type) )
#AFN <- subset_taxa(AFN, Phylum != "" & Genus != "" & !is.na(Genus))
#meli <- psmelt(AFN) %>% filter(Abundance > 1000) %>% select(Abundance, ecoregion.III,env.broad.scale,animal.type, Diet, Phylum, Genus)
#alluvial_plot <- alluvial_plot(data = meli, axes = c("animal.type", "ecoregion.III", "Diet"), abundance_threshold = 1000, silent = TRUE)
#alluvial_plot
