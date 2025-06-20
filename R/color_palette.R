#' @title Original, Extended, and Nature-Inspired Color Palette Sequence
#' @return A named list of character vectors, each containing color hex codes for different palettes.
#' @description This object provides multiple color palettes for scientific, exploratory, and publication-grade graphics:
#'  - MG: The original palette.
#'  - extended_palette: MG combined with rainbow colors.
#'  - light_MG: A chic pastel palette.
#'  - MG_Awesome: A vibrant and unique 50-color palette for visualization.
#'  - cool_MG: A sophisticated, modern 50-color palette with oceanic, earthy, and high-contrast tones.
#'  - mix_MG: A randomly mixed unique palette.
#'  - vivid_MG: A vivid and bright 50-color palette.
#'  - Mar_palette: A carefully crafted 40-color nature-inspired palette suitable for high-quality publications.
#'
#' @details These palettes can be directly used with ggplot2 or other visualization systems supporting manual color scales.
#'
#' @format A list with eight elements:
#' \describe{
#'   \item{MG}{A character vector of the original color codes.}
#'   \item{extended_palette}{A character vector of extended color codes, combining MG with the rainbow palette.}
#'   \item{light_MG}{A character vector of chic pastel color codes.}
#'   \item{MG_Awesome}{A vibrant 50-color palette designed for high-impact visualizations.}
#'   \item{cool_MG}{A sophisticated and modern 50-color palette with oceanic, earthy, and high-contrast tones.}
#'   \item{mix_MG}{A fully combined and randomly mixed unique palette including Mar_palette.}
#'   \item{vivid_MG}{A vivid and bright 50-color palette.}
#'   \item{Mar_palette}{A 40-color palette inspired by nature, suitable for academic and professional publications.}
#' }
#'
#' @examples
#' # Example using the Mar_palette
#' ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, color = factor(cyl))) +
#'   ggplot2::geom_point(size = 4) +
#'   ggplot2::scale_color_manual(values = color_palette$Mar_palette) +
#'   ggplot2::theme_minimal()
#'
#' @name color_palette
#' @importFrom grDevices rainbow
#' @export
color_palette <- local({
  # Define the original MG color palette
  MG <- c(
    "#FF7F00", "#E41A1C", "#BB650B", "firebrick4", "#FC4E07", "gold", "#E7B800", "#FF9933",
    "#CC9900", "#AE2012", "#BB3E03", "#CA6702", "#EE9B00", "#9B2226", "brown", "darkgoldenrod",
    "#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#2e4057", "#377EB8", "#66C2A5", "#0e669b",
    "#00798c", "dodgerblue4", "steelblue2", "#00AFBB", "#0099CC", "#0000CC", "#17b5b4", "#AED1D6",
    "#82cfd0", "#b2e0e4", "lightskyblue4", "#0e668b", "#a3c4dc", "lightskyblue1", "aliceblue",
    "olivedrab3", "#4DAF4A", "#336633", "#66a61e", "chartreuse1", "#A6D854", "#909800", "#00FF00",
    "green", "darkgreen", "#386641", "#1B4332", "#2D6A4F", "#40916C", "#52B788", "#74C69D",
    "#984EA3", "#E78AC3", "#7570b3", "#9183E6", "#00AD9A", "#CC6686", "#FF3399", "#990033",
    "#6A0572", "#AB83A1", "#553C9A", "#3F37C9", "#8338EC", "#DA33FF", "purple", "darkred",
    "#663300", "#A37F6F", "#e1deda", "burlywood", "papayawhip", "wheat4", "cornsilk3", "khaki2",
    "beige", "gray60", "gray80", "gray96", "cadetblue4", "honeydew2", "mintcream", "#8d96a3",
    "lavender", "lavenderblush2", "mistyrose3", "black", "navy", "yellow"
  )

  # Extended palette using the rainbow function
  extended_palette <- c(MG, grDevices::rainbow(50))

  # Chic Pastel Palette (light_MG)
  light_MG <- c(
    "#F1E0C5", "#D2A5A1", "#B2C3A8", "#B8B1A3", "#A0869A", "#C4AB86", "#4F4A4A", "#FFD2A5",
    "#F6E7D8", "#BFC4C9", "#8097A4", "#A3AC9A", "#D0D3C3", "#F1CAB8", "#D49F9B", "#C2AB99",
    "#9A8475", "#A9B7B1", "#C7D1CC", "#6F7367", "#F4C8BD", "#D49C8E", "#C8ABAD", "#8D8885",
    "#E4D9D2", "#A79A91", "#B1B6A1", "#D2CFC1", "#F1E5E5", "#BCB1A3", "#D3A3AB", "#E0C7B6",
    "#998C8B", "#A9A59B", "#B8B3A3", "#8D897C", "#C3BEB5", "#A5938E", "#E6C7C1", "#C1B1A9",
    "#E8D3D1", "#D0AFA5", "#B5AAA5", "#E0DED9", "#B0A59D", "#C7B8AE", "#9F9894", "#D5D0CB",
    "#CCC6C1", "#E2DAD7"
  )

  # Vibrant and Unique Palette (MG_Awesome)
  MG_Awesome <- c(
    "#FF5733", "#33FF57", "#3357FF", "#FF33A8", "#A833FF",
    "#33FFF5", "#F5FF33", "#FF8C33", "#33FF8C", "#8C33FF",
    "#FF336E", "#6EFF33", "#336EFF", "#FFDA33", "#DA33FF",
    "#33FFD1", "#D133FF", "#FFD133", "#33B8FF", "#B833FF",
    "#E63946", "#F4A261", "#2A9D8F", "#264653", "#E76F51",
    "#1B998B", "#ED217C", "#2D3047", "#FF9B71", "#FFD700",
    "#17BEBB", "#EF476F", "#FFD166", "#06D6A0", "#118AB2",
    "#073B4C", "#F94144", "#F3722C", "#F8961E", "#43AA8B",
    "#277DA1", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
    "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548"
  )

  # Sophisticated Cool Colors Palette
  cool_MG <- c(
    "#001219", "#CA6702", "#780000", "#E9C46A", "#386641", "#669BBC",
    "#005F73", "#BB3E03", "#C1121F", "#F4A261", "#1B4332", "#003049",
    "#0A9396", "#AE2012", "#FDF0D5", "#E76F51", "#2D6A4F", "#8338EC",
    "#94D2BD", "#9B2226", "#553C9A", "#6A0572", "#40916C", "#FF006E",
    "#E9D8A6", "#582F0E", "#3F37C9", "#AB83A1", "#52B788", "#FB5607",
    "#EE9B00", "#7F4F24", "#4895EF", "#FCA311", "#74C69D", "#FFBE0B",
    "#264653", "#936639", "#14213D", "#E5E5E5", "#8D99AE", "#3A86FF",
    "#2A9D8F", "#A68A64", "#FFBA08", "#EDF2F4", "#D00000", "#BC4749",
    "#6A994E", "#B6AD90", "#8D99AE", "#D00000", "#FF006E", "#780000"
  )

  # vivid colors

  vivid_MG <- c(
    "#FF0000", "#FF4500", "#FF8C00", "#FFD700", "#FFFF00", "#ADFF2F", "#7CFC00", "#32CD32",
    "#00FF00", "#00FA9A", "#00FFFF", "#00CED1", "#1E90FF", "#4169E1", "#0000FF", "#8A2BE2",
    "#9400D3", "#FF00FF", "#EE82EE", "#DA70D6", "#FFC0CB", "#FF69B4", "#FF1493", "#DC143C",
    "#B22222", "#FA8072", "#E9967A", "#FFA07A", "#FF4500", "#FF6347", "#D2691E", "#8B4513",
    "#FFDAB9", "#FFDEAD", "#F4A460", "#D2B48C", "#BDB76B", "#ADFF2F", "#9ACD32", "#6B8E23",
    "#556B2F", "#008B8B", "#20B2AA", "#48D1CC", "#40E0D0", "#5F9EA0", "#4682B4", "#6495ED",
    "#7B68EE", "#8B0000"
  )

  # MarMar
  Mar_palette <- c(
    "#1b4332", "#2d6a4f", "#386641", "#40916c", "#52b788",
    "#74c69d", "#90be6d", "#b5c99a", "#ccd5ae", "#e9d8a6",
    "#fefae0", "#faedcd", "#f4a261", "#e76f51", "#bc6c25",
    "#dda15e", "#d4a373", "#c19a6b", "#a98467", "#8c6f58",
    "#6c584c", "#5f4b3a", "#b08968", "#e0ac9d", "#c4a69f",
    "#b99382", "#997b66", "#d7b29d", "#ecd9c6", "#ffe8d6",
    "#6a994e", "#a7c957", "#f9bc60", "#e36414", "#9a031e",
    "#5f0f40", "#0a9396", "#005f73", "#ee9b00", "#3f37c9"
  )

  mix_MG <- unique(sample(c(cool_MG, light_MG, MG, MG_Awesome, vivid_MG, Mar_palette)))

  list(
    Mar_palette = Mar_palette,
    MG = MG,
    extended_palette = extended_palette,
    light_MG = light_MG,
    MG_Awesome = MG_Awesome,
    cool_MG = cool_MG,
    vivid_MG = vivid_MG,
    mix_MG = mix_MG
  )
})
# 
