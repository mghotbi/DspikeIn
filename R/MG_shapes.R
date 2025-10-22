#' @title Predefined Shape Vector for Plot Styling
#' @description A numeric vector of point shapes (0-25) for use in ggplot2.
#' @return A numeric vector of shape codes.
#'
#' @examples
#' # Example usage of MG_shapes in ggplot2
#' ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, shape = factor(cyl))) +
#'   ggplot2::geom_point(size = 4) +
#'   ggplot2::scale_shape_manual(values = MG_shapes[seq_len(3)]) +
#'   my_custom_theme()
#'
#' @export
MG_shapes <- c(
  21, 24, 22, 23, 25, 0, 3, 4, 8, 14, 13, 12, 11, 9, 10, 6, 7, 15,
  16, 17, 18, 19, 20, 2, 5, 1
)
