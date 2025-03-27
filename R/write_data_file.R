#' Write a data file for the AgeingError package
#'
#' The resulting file has the format required by the [load_data()] function.
#' In the future, this step could be bypassed by creating the list output
#' from [load_data()] directly. Only one data set is supported (NDataSet = 1).
#' This function is based on a subset of the RunFn() function in the older ADMB
#' version of this package.
#'
#' @param dat Dataframe or tibble with columns for each reader and rows for
#'   each age reading combination. This could either have a count column or
#'   not. If not, [tally_repeats()] will be called to add a count column.
#'   Order your reader/lab columns such that similar
#'   readers/labs are located next to one another because columns to the right
#'   can mirror columns to their immediate left in terms of parameter
#'   estimates.
#' @param dir Directory where the data file will be saved.
#' @param file_name Name of the data file.
#' @param minage An integer, specifying the minimum possible "true" age.
#' @param maxage An integer, specifying the maximum possible "true" age.
#' @param refage An arbitrarily chosen age from which "true" age-composition
#'   fixed-effects are calculated as an offset. This has no effect on the
#'   answer but could potentially effect estimation speed. By default this will
#'   be set to the maxage / 4.
#' @param minusage The minimum age for which an age-specific age-composition is
#'   estimated. Ages below `minusage` have "true" proportion-at-age
#'   (\eqn{P_{a}}) estimated as
#'   \deqn{P_a = P_{minusage}*exp^{(\beta*(minusage - a))}},
#'   where beta is an estimated log-linear trend in the "true"
#'   proportion-at-age. If `minusage` = `minage`, beta is not estimated.
#' @param plusage Identical to `minusage` except defining the age above with
#' age-specific age composition is not estimated.
#' @return Invisibly returns the path to the data file (`file.path(dir, file_name)`).
#' @author Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt
#' @export
#' @seealso [write_files()], [write_specs_file()], [load_data()], [tally_repeats()]
#' @examples
#' data_test <- data.frame(
#'   reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
#'   reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
#'   reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
#' )
#' data_file <- write_data_file(data_test, dir = tempdir(), file_name = "test.dat")
#'
write_data_file <- function(
    dat,
    dir = getwd(),
    file_name = "data.dat",
    minage = 0,
    maxage = NULL,
    refage = NULL,
    minusage = NULL,
    plusage = NULL) {
  # check inputs
  if (!is.data.frame(dat)) {
    cli::cli_abort("Input 'dat' must be a data frame or tibble")
  }
  if (!dir.exists(dir)) {
    cli::cli_alert_info("Directory does not exist; creating it")
    dir.create(dir)
  }
  if (!is.character(file_name)) {
    cli::cli_abort("Input 'file_name' must be a character string")
  }

  # add count column if not already present
  if (!"count" %in% colnames(dat)) {
    cli::cli_alert_info(
      "Input 'dat' doesn't contain a column called 'count'; adding one via tally_repeats()"
    )
    dat <- tally_repeats(dat)
  }

  # replace any NA values with -999 otherwise the tally gets messed up
  if (any(is.na((dat)))) {
    test_data[is.na(dat)] <- -999
  }

  minobs <- min(abs(dat[, -1])) # abs is to avoid counting the -999 values
  maxobs <- max(dat[, -1])
  nreaders <- ncol(dat) - 1
  cli::cli_alert_info("Range of observed ages in the data: {minobs} - {maxobs}")
  cli::cli_alert_info("Number of readers: {nreaders}")

  # fill in any missing default values
  if (is.null(maxage)) {
    maxage <- ceiling(1.2 * maxobs / 5) * 5
    cli::cli_alert_info(
      "Max age not specified; using {maxage} which is the multiple of 5 which is >120% of the observed maximum"
    )
  }
  if (is.null(minusage)) {
    minusage <- minobs
    cli::cli_alert_info(
      "Minus group set to the minimum observed age {minusage}"
    )
  }
  if (is.null(plusage)) {
    plusage <- maxobs
    cli::cli_alert_info("Plus group set to the maximum observed age {plusage}")
  }
  if (is.null(refage)) {
    refage <- floor(median(c(minusage, plusage)))
    cli::cli_alert_info(
      "Reference age not specified; using {refage} = floor(median(c(minusage, plusage)))"
    )
  }

  # create the header
  header <- c(
    "Range_of_ages",
    paste(minage, maxage),
    "",
    "Data_set_1",
    paste(nrow(dat), "# data points"),
    paste(nreaders, "# readers"),
    paste(
      minusage,
      plusage,
      refage,
      "# minus group; plus group; reference age"
    ),
    paste0("   ", 1:nreaders, collapse = " ")
  )

  # write file
  cli::cli_alert_info("Writing data file to {file.path(dir, file_name)}")
  # write the header
  writeLines(
    text = header,
    con = file.path(dir, file_name)
  )

  # write data
  utils::write.table(
    dat,
    file = file.path(dir, file_name),
    append = TRUE,
    row.names = FALSE,
    col.names = FALSE
  )

  return(invisible(file.path(dir, file_name)))
}
