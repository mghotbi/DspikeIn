#' @title Create a Directory and Optionally Set as Working Directory
#'
#' @description This function checks if a specified directory exists and creates it if it doesn't.
#' Optionally, it can also set the newly created or existing directory as the working directory.
#'
#' @param directory_path A character string specifying the path of the directory to create.
#' @param set_working_dir A logical value indicating whether to set the directory as the working directory. Default is FALSE.
#' @return NULL. The function prints messages indicating whether the directory was created or already exists, and if the working directory was set.
#' @examples
#' if (interactive()) {
#'   # Save the current working directory
#'   old_wd <- getwd()
#'
#'   # Use a temporary directory for safe example use
#'   tmp_dir <- file.path(tempdir(), "example_new_dir")
#'
#'   # Create the directory and set it as the working directory
#'   create_directory(tmp_dir, set_working_dir = TRUE)
#'
#'   # Do something inside tmp_dir...
#'
#'   # Restore the original working directory
#'   setwd(old_wd)
#'
#'   # Remove the created directory
#'   unlink(tmp_dir, recursive = TRUE, force = TRUE)
#' }
#' @export
create_directory <- function(directory_path, set_working_dir = FALSE) {
  if (!base::dir.exists(directory_path)) {
    base::dir.create(directory_path)
    message("Directory created successfully.")
  } else {
    message("Directory already exists.")
  }

  if (set_working_dir) {
    setwd(directory_path)
    message("Working directory set to: ", directory_path)
  }
}


# Example usage:
# create_directory("new", set_working_dir = TRUE)
