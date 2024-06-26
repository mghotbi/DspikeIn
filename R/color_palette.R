# MG color_palettes

# Original color palette sequence
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

# Extend the color palette sequence
extended_palette <- c(MG, rainbow(50))

# Export the color palette
color_palette <- list(MG = MG, extended_palette = extended_palette)

#' Original and Extended Color Palette Sequence
#'
#' This object provides an original color palette sequence (MG) and an extended palette using the rainbow function.
#' The extended palette can be used for plotting in ggplot2 or other plotting systems.
#'
#' @format A list of character vectors of color codes.
#' @examples
#' # Use the extended palette for your plot
#' # ggplot2 example:
#' # ggplot(data, aes(x, y, color = group)) +
#' #   geom_point() +
#' #   scale_color_manual(values = extended_palette)
#' @name color_palette
#' @export
NULL
