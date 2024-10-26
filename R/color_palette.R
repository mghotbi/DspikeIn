#' Original and Extended Color Palette Sequence
#'
#' This object provides an original color palette sequence (MG), an extended palette using the rainbow function,
#' and an additional chic palette named light_MG.
#' These palettes can be used for plotting in ggplot2 or other plotting systems.
#'
#' @format A list with three elements: 
#' \describe{
#'   \item{MG}{A character vector of original color codes.}
#'   \item{extended_palette}{A character vector of extended color codes, combining the original palette with the rainbow palette.}
#'   \item{light_MG}{A character vector of chic palette color codes.}
#' }
#' @examples
#' # Use the light_MG palette for your plot
#' # ggplot2 example:
#' # ggplot(data, aes(x, y, color = group)) +
#' #   geom_point() +
#' #   scale_color_manual(values = color_palette$light_MG)
#' @name color_palette
#' @importFrom grDevices rainbow
#' @export
color_palette <- local({
  MG <- c("#FFFF33", "#FF7F00", "#E41A1C", "firebrick4", "#2e4057", "#984EA3", "#377EB8", 
          "olivedrab3", "#4DAF4A", "#336633", "grey80", "#BB650B", "gold", "#559999", 
          "#7570b3", "#E78AC3", "#A6D854", "#66a61e", "#e6ab02", "#a6761d", "#663300", 
          "#66C2A5", "#0e669b", "#00798c", "dodgerblue4", "steelblue2", "#00AFBB", 
          "#E7B800", "#FC4E07", "lightskyblue4", "green", "red", "#FFF000", "#0099CC", 
          "#FF9933", "#CC9900", "chartreuse1", "#FF3399", "#00FFFF", "#0000CC", "#A37F6F", 
          "#9183E6", "#00AD9A", "#990033", "#909800", "#00FF00", "#17b5b4", "#AED1D6", 
          "#b1010c", "firebrick2", "blue", "navy", "yellow", "brown", "black", "purple", 
          "darkred", "darkgreen", "#82cfd0", "#b2e0e4", "honeydew3", "#8d96a3", "lavender", 
          "#CC6686", "lavenderblush2", "mistyrose3", "#e1deda", "darkgoldenrod", "burlywood", 
          "papayawhip", "wheat4", "cornsilk3", "khaki2", "beige", "gray60", "gray80", 
          "gray96", "cadetblue4", "honeydew2", "mintcream", "#0e668b", "#a3c4dc", 
          "lightskyblue1", "aliceblue")
  
  # Create the extended palette using the rainbow function
  extended_palette <- suppressMessages(grDevices::rainbow(50))
  
  # Define the chic_palette as light_MG
  chic_palette <- c(
    "#F1E0C5", "#D2A5A1", "#B2C3A8", "#B8B1A3", "#A0869A", "#C4AB86", "#4F4A4A", "#FFD2A5", 
    "#F6E7D8", "#BFC4C9", "#8097A4", "#A3AC9A", "#D0D3C3", "#F1CAB8", "#D49F9B", "#C2AB99", 
    "#9A8475", "#A9B7B1", "#C7D1CC", "#6F7367", "#F4C8BD", "#D49C8E", "#C8ABAD", "#8D8885", 
    "#E4D9D2", "#A79A91", "#B1B6A1", "#D2CFC1", "#F1E5E5", "#BCB1A3", "#D3A3AB", "#E0C7B6", 
    "#998C8B", "#A9A59B", "#B8B3A3", "#8D897C", "#C3BEB5", "#A5938E", "#E6C7C1", "#C1B1A9", 
    "#E8D3D1", "#D0AFA5", "#B5AAA5", "#E0DED9", "#B0A59D", "#C7B8AE", "#9F9894", "#D5D0CB", 
    "#CCC6C1", "#E2DAD7"
  )
  
  # Return the list with all three palettes
  list(
    MG = MG,
    extended_palette = c(MG, extended_palette),
    light_MG = chic_palette  # Add the chic_palette as light_MG
  )
})

# Example usage:
# print(color_palette$MG)  # Prints the original MG palette
# print(color_palette$extended_palette)  # Prints the extended palette
# print(color_palette$light_MG) # Prints the light MG
