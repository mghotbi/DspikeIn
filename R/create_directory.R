#' @title Create a Directory and Optionally Set as Working Directory
#'
#' @description This function checks if a specified directory exists and creates it if it doesn't.
#' Optionally, it can also set the newly created or existing directory as the working directory.
#'
#' @param directory_path A character string specifying the path of the directory to create.
#' @param set_working_dir A logical value indicating whether to set the directory as the working directory. Default is FALSE.
#' @return NULL. The function prints messages indicating whether the directory was created or already exists, and if the working directory was set.
#' @examples
#' # Create a new directory and set it as the working directory
#' create_directory("new", set_working_dir = TRUE)
#'
#' # Clean up afterwards
#' setwd(tempdir())
#' unlink("new", recursive = TRUE)
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
