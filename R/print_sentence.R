#' Print a Sentence Multiple Times
#'
#' This function prints a specified sentence multiple times (default is 10 times).
#'
#' @param sentence A character string specifying the sentence to be printed.
#' @return NULL. The function prints the sentence to the console.
#' @examples
#' # Example usage:
#' # Define a sentence
#' sentence <- " ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯ 
#' # ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯ "
#'
#' # Print the sentence multiple times
#' print_sentence(sentence)
#' @export
print_sentence <- function(sentence) {
  for (i in 1:10) { 
    cat(sentence, "\n")
  }
}

# Example usage:
# Define a sentence
#sentence <- " ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯  
#¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯ ¯\\_(ツ)_/¯  ¯\\_(ツ)_/¯ "

# Print the sentence multiple times
#print_sentence(sentence)
