# File: zzz.R

#' @title Package Startup Message
#' @description Displays a message when the DspikeIn package is loaded.
#' @param libname The library name (automatically passed by R).
#' @param pkgname The package name (automatically passed by R).
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  msg <- paste0(
    "\n",
    "Thank you for using **DspikeIn**! \U0001F389\n",
    "For support,\U0001F4E7 contact: Mitra Ghotbi (mitra.ghotbi@gmail.com)\n",
    " GitHub Repository: https://github.com/mghotbi/DspikeIn \n",
    "\U0001F41B Found a bug or have a suggestion? Open an issue on GitHub: \n",
    "   https://github.com/mghotbi/DspikeIn/issues \n",
    "\n"
  )
  packageStartupMessage(msg)
}
