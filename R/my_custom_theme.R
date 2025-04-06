#' @title Custom ggplot2 Theme with Consistent Aesthetics
#' @description Creates a custom ggplot2 theme with consistent styling options,
#' including background color, font size, and axis line formatting.
#'
#' @param base_size Numeric. Base font size for text elements (default: 12).
#' @param font_family Character. Font family to use for text elements (default: "sans").
#' @param bg_color Character. Background color of the plot (default: "white").
#' @return A ggplot2 theme object.
#'
#' @details
#' This function applies a custom theme to ggplot2 plots. It removes unnecessary
#' grid lines, ensures a clean background, and standardizes text styling.
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point(size = 3) +
#'   my_custom_theme()
#' print(p)
#'
#' @importFrom ggplot2 theme_bw theme element_blank element_rect element_line element_text unit
#' @export
my_custom_theme <- function(base_size = 12, font_family = "sans", bg_color = "white") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = bg_color),
      plot.background = ggplot2::element_rect(fill = bg_color, color = "#e1deda"),
      panel.border = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(colour = "black", size = 0.6),
      axis.line.y = ggplot2::element_line(colour = "black", size = 0.6),
      axis.ticks = ggplot2::element_line(colour = "black", size = 0.35),
      legend.title = ggplot2::element_text(size = base_size, family = font_family),
      legend.text = ggplot2::element_text(size = base_size - 1, color = "black", family = font_family, face = "plain"),
      legend.key.size = ggplot2::unit(0.7, "cm"),
      axis.title.x = ggplot2::element_text(size = base_size, color = "black", family = font_family, face = "bold"),
      axis.title.y = ggplot2::element_text(size = base_size, color = "black", family = font_family, face = "bold"),
      axis.text.x = ggplot2::element_text(size = base_size, color = "black", family = font_family, face = "bold"),
      axis.text.y = ggplot2::element_text(size = base_size, color = "black", family = font_family, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", family = font_family),
      plot.subtitle = ggplot2::element_text(size = base_size, family = font_family),
      strip.background = ggplot2::element_rect(fill = "gray98", color = NA),
      strip.text = ggplot2::element_text(size = base_size + 1, face = "bold", family = font_family)
    )
}
