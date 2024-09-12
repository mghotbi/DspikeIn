#' Predefined shape vector for consistent plot styling
#'
#' This vector provides predefined shapes for consistent plot styling in ggplot2 
#'
#' @usage MG_shape
#' @format A numeric vector.
#' @examples
#' # Example usage of MG_shape in a base R plot
#' plot(1:10, pch = MG_shape[1:10]) # Use first 10 shapes from MG_shape
#'
#' # Example usage in ggplot2
#' # ggplot(mtcars, aes(x = wt, y = mpg, shape = factor(cyl))) + 
#' # geom_point(size = 4) + scale_shape_manual(values = MG_shape[1:3])
#' @export
MG_shape <- c(19, 3, 1, 2, 9, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 21, 23, 20, 22)
