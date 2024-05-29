#' Determine the number of data sets in a data file
#' @param file A file path to a data file.
#' @author Kelli F. Johnson
#' @return An integer giving the number of data sets in the file.
determine_n_sets <- function(file) {
  the_lines <- readLines(file)
  # Strip down to numbers and spaces
  the_numbers <- gsub("\\p+|[[:alpha:][:punct:]]+", "", the_lines)
  the_end <- gsub("([0-9])\\s+$", "\\1", the_numbers)
  single_numbers <- grep("^[0-9]+$", the_end)
  number_of_sets <- length(single_numbers) / 3
  stopifnot(number_of_sets %% 1 == 0)
  return(as.integer(number_of_sets))
}
