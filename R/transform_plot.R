#' Generate and Save Boxplots with Statistical Tests
#'
#' This function generates boxplots for multiple y variables and performs statistical tests (Kruskal-Wallis or ANOVA) to compare groups. It also generates individual boxplots for Y variables and performs statistical tests (Kruskal-Wallis or Wilcoxon) to compare groups, and saves the plots as PNG and PDF files.
#'
#' @param X A numeric vector of data values.
#' @param Y A factor or character vector of grouping variables.
#' @param data A data frame containing the data to plot.
#' @param x_var A character string specifying the column name for the x variable.
#' @param y_vars A character vector specifying the column names for the y variables.
#' @param methods_var A character string specifying the column name for the grouping variable.
#' @param MG A character vector specifying the colors for the boxplots. Default is MG.
#' @param adjustment A character string specifying the method for p-value adjustment. Default is "holm".
#' @param output_prefix A character string specifying the prefix for the output file names. Default is "plot".
#' @param width A numeric value specifying the width of the output plot. Default is 15.
#' @param height A numeric value specifying the height of the output plot. Default is 13.
#' @param stat_test A character string specifying the statistical test to use ("kruskal.test" or "anova"). Default is "kruskal.test".
#' @param main A character string specifying the main title of the plot. Default is NULL.
#' @param xlab A character string specifying the x-axis label. Default is NULL.
#' @param ylab A character string specifying the y-axis label. Default is NULL.
#' @param bcol A character string specifying the box color. Default is "bisque".
#' @param p.adj A character string specifying the method for p-value adjustment. Default is "none".
#' @param cexy A numeric value specifying the text size for axis labels and titles. Default is 1.5.
#' @param varwidth A logical value specifying whether the boxes should have variable widths. Default is TRUE.
#' @param las An integer specifying the orientation of axis labels. Default is 1.
#' @param paired A logical value specifying whether the data are paired. Default is FALSE.
#' @return A list containing the ggplot2 boxplot objects and the comparison results, or the comparison results and p-value for individual boxplots.
#' @examples
#' # Example usage for transform_plot:
#' # y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")
#' # transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", MG = MG, stat_test = "kruskal.test")
#' # transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", MG = MG, stat_test = "anova")
#' #
#' # Example usage for individual boxplots:
#' # X <- rnorm(100)
#' # Y <- rep(c("A", "B"), each = 50)
#' # transform_plot(X = X, Y = Y)
#' @export
transform_plot <- function(
    data = NULL,
    x_var = NULL,
    y_vars = NULL,
    methods_var = NULL,
    MG = MG,
    adjustment = "holm",
    output_prefix = "plot",
    width = 15,
    height = 13,
    stat_test = "kruskal.test",
    X = NULL,
    Y = NULL,
    main = NULL,
    xlab = NULL,
    ylab = NULL,
    bcol = "bisque",
    p.adj = "none",
    cexy = 1.5,
    varwidth = TRUE,
    las = 1,
    paired = FALSE
) {
  # Load necessary libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required but not installed.")
  }
  
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  
  if (!is.null(data) && !is.null(x_var) && !is.null(y_vars) && !is.null(methods_var) && !is.null(MG)) {
    # Suffix for file names based on stat_test
    test_suffix <- ifelse(stat_test == "kruskal.test", "Kruskal", "ANOVA")
    
    # Function to create individual boxplots
    create_boxplot <- function(y_var) {
      ggplot(data, aes(x = .data[[methods_var]], y = .data[[y_var]], fill = .data[[methods_var]])) +
        ggplot2::geom_boxplot() +
        scale_fill_manual(values = MG) +
        labs(x = "", y = y_var, title = y_var) +
        theme_minimal() +
        theme(
          axis.text = element_text(size = 20, angle = 20, hjust = 1),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22)
        ) +
        stat_compare_means(method = stat_test, label = "p.signif", size = 8)
    }
    
    # Create and save plots
    for (y_var in y_vars) {
      plot <- create_boxplot(y_var)
      print(plot)
      
      png_filename <- paste0(output_prefix, "_", y_var, "_", test_suffix, ".png")
      pdf_filename <- paste0(output_prefix, "_", y_var, "_", test_suffix, ".pdf")
      
      # Save PDF
      ggsave(pdf_filename, plot = plot, width = width, height = height, units = "in")
      # Save PNG with explicit DPI to ensure text size is consistent
      ggsave(png_filename, plot = plot, width = width, height = height, units = "in", dpi = 500)
      
      cat("Plots saved as:", png_filename, "and", pdf_filename, "\n")
    }
  }
  
  if (!is.null(X) && !is.null(Y)) {
    aa <- levels(as.factor(Y))
    an <- as.character(c(1:length(aa)))
    tt1 <- matrix(nrow = length(aa), ncol = 7)
    
    for (i in 1:length(aa)) {
      temp <- X[Y == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- median(temp, na.rm = TRUE)
      tt1[i, 7] <- length(temp)
    }
    
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
    
    boxplot(
      X ~ Y,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth,
      cex.main = cexy + 0.5
    )
    
    Yn <- factor(Y, labels = an)
    comp <- kruskal(X, Yn, p.adj = p.adj)
    sig <- "ns"
    
    if (paired == TRUE & length(aa) == 2) {
      coms <- wilcox.test(X ~ Yn, paired = TRUE)
      pp <- coms$p.value
    } else {
      pp <- comp$statistics$p.chisq
    }
    
    if (pp <= 0.1) sig <- "."
    if (pp <= 0.05) sig <- "*"
    if (pp <= 0.01) sig <- "**"
    if (pp <= 0.001) sig <- "***"
    
    gror <- comp$groups[order(rownames(comp$groups)), ]
    tt1$rank <- gror$X
    tt1$group <- gror$groups
    mtext(
      sig,
      side = 3,
      line = 1,
      adj = 0,
      cex = 3,
      font = 1
    )
    if (pp <= 0.1)
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 1,
        cex = 2,
        font = 4
      )
    
    return(list(comparison = tt1, p.value = pp))
  }
}

# Example y_vars
#y_vars <- c("Spike.percentage", "Total.reads", "Spike.reads")

# Ensure the columns are numeric
#methods <- methods %>%
 # dplyr::mutate(
  #  Total.reads = as.numeric(Total.reads),
  #  Spike.reads = as.numeric(Spike.reads),
   # Spike.percentage = as.numeric(Spike.percentage)  )

#print(sapply(methods[, c("Total.reads", "Spike.reads", "Spike.percentage")], class))

# Remove rows with NA values
#methods <- methods %>%
  #dplyr::filter(
  #  !is.na(Total.reads),
   # !is.na(Spike.reads),
   # !is.na(Spike.percentage)  )

# Scale the specified columns
#scaled <- methods %>%
#  dplyr::mutate_at(
 #   c("Total.reads", "Spike.reads", "Spike.percentage"),    ~ scale(.) %>% as.vector  )

# Perform Kruskal-Wallis test
#transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", colors = MG, stat_test = "kruskal.test")

# Perform one-way ANOVA
#transform_plot(data = scaled, x_var = "Methods", y_vars = y_vars, methods_var = "Methods", colors = MG, stat_test = "anova")

# Perform individual boxplots with Kruskal-Wallis test
#X <- rnorm(100)
#Y <- rep(c("A", "B"), each = 50)
#transform_plot(X = X, Y = Y)
