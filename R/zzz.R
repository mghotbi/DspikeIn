# File: zzz.R

#' Startup message for DspikeIn package
#'
#' This function displays a startup message when the DspikeIn package is loaded.
#'
#' @param libname The library name.
#' @param pkgname The package name.
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage('Thank you for using the DspikeIn package.')
  packageStartupMessage('For support, please contact Mitra Ghotbi at mitra.ghotbi@gmail.com.')
}
