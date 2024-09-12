#' Custom GGplot2 Theme and Color/Shape Vectors
#'
#' This function generates a custom ggplot2 theme with specific styling options.
#' Additionally, it provides predefined color and shape vectors for consistent plot styling.
#'
#' @return A ggplot2 theme object.
#' @importFrom ggplot2 theme_bw theme element_blank element_rect element_line element_text unit
#' @examples
#' # Apply the custom theme to a ggplot
#' suppressMessages(library(ggplot2))
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) + 
#'   geom_point(size = 3) + 
#'   my_custom_theme()
#' print(p)
#' @export
my_custom_theme <- function() {
  ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background = ggplot2::element_rect(fill = 'white', color = "#e1deda"),
      panel.border = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = 'black', size = 0.6),
      axis.line.y = ggplot2::element_line(colour = 'black', size = 0.6),
      axis.ticks = ggplot2::element_line(colour = 'black', size = 0.35),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11, color = "black", face = "plain"), 
      legend.key.size = ggplot2::unit(0.7, 'cm'),
      axis.title.x = ggplot2::element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), 
      axis.title.y = ggplot2::element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), 
      axis.text.x = ggplot2::element_text(family = "Times New Roman", size = 12, angle = 0, color = "black", face = "bold"), 
      axis.text.y = ggplot2::element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
      plot.title = ggplot2::element_text(color = "black", size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      strip.background = ggplot2::element_rect(fill = "gray98", color = NA),  # Set facet background color
      strip.text = ggplot2::element_text(size = 13, face = "bold")  # Increase facet text size
    )
}

#' Predefined color and shape vectors for consistent plot styling
#' @export
MG <- c("#FFFF33", "#FF7F00", "#E41A1C", "firebrick4", "#2e4057", "#984EA3", "#377EB8", "olivedrab3", 
        "#4DAF4A", "#336633", "grey80", "#BB650B", "gold", "#559999", "#7570b3", "#E78AC3", "#A6D854", 
        "#66a61e", "#e6ab02", "#a6761d", "#663300", "#66C2A5", "#0e669b", "#00798c", "dodgerblue4", 
        "steelblue2", "#00AFBB", "#E7B800", "#FC4E07", "lightskyblue4", "green", "red", "#FFF000", 
        "#0099CC", "#FF9933", "#CC9900", "chartreuse1", "#FF3399", "#00FFFF", "#0000CC", "#A37F6F", 
        "#9183E6", "#00AD9A", "#990033", "#909800", "#00FF00", "#17b5b4", "#AED1D6", "#b1010c", 
        "firebrick2", "blue", "navy", "yellow", "brown", "black", "purple", "darkred", "darkgreen", 
        "#82cfd0", "#b2e0e4", "honeydew3", "#8d96a3", "lavender", "#CC6686", "lavenderblush2", 
        "mistyrose3", "#e1deda", "darkgoldenrod", "burlywood", "papayawhip", "wheat4", "cornsilk3", 
        "khaki2", "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream", 
        "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue")

#' @export
MG_shape <- c(19, 3, 1, 2, 9, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 23, 20, 22)
