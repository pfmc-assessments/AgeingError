#' Tally repeated age reading combinations in a new count column
#'
#' AgeingError requires a table of data with columns for each reader and
#' rows for each age reading combination, where a `count` column indicates
#' how many times the age reading combination is repeated. This function
#' tallies the repeated values and adds an initial `count` column with the
#' tally.
#'
#' @param dat Input dataframe or tibble with columns for each reader and rows
#' for each specimen
#' @export
#' @author Ian G. Taylor

tally_repeats <- function(dat) {
  # add dummy column with value 1 for all rows if not already present
  if (!"count" %in% colnames(dat)) {
    dat$count <- 1
  }

  # replace any NA values with -999 otherwise the tally gets messed up
  if (any(is.na((dat)))) {
    dat[is.na(dat)] <- -999
  }

  # aggregate duplicates to change count
  dat2 <- dat |>
    dplyr::group_by_all() |>
    dplyr::summarize(count = sum(count)) |>
    dplyr::relocate(count, 1) # put count in the first column

  # messages about rows of input and output data
  cli::cli_alert_info("Total observations: {nrow(dat)}")
  cli::cli_alert_info("Aggregated unique combinations: {nrow(dat2)}")

  return(dat2)
}
