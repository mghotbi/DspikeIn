#' @title Load FASTA Sequence Data
#'
#' @description This function loads a FASTA file stored in `inst/extdata/` of the package.
#'
#' @param filename Name of the FASTA file (default: `"combined.fasta"`).
#' @return A `DNAStringSet` object containing the sequences.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("DspikeIn", quietly = TRUE)) {
#'   # Load FASTA sequences from a file
#'   sequences <- load_fasta("combined.fasta")
#' }
#' }
#' @importFrom Biostrings readDNAStringSet
#' @export
load_fasta <- function(filename = "combined.fasta") {
  file_path <- system.file("extdata", filename, package = "DspikeIn")

  if (file_path == "") {
    stop("Error: File not found in the installed package's inst/extdata/")
  }

  Biostrings::readDNAStringSet(file_path)
}
