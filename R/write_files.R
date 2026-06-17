#' Write the input files (data and specifications) for the AgeingError package
#'
#' @inheritParams write_data_file
#' @inheritParams write_specs_file
#' @return Invisibly returns the path to the data file (`file.path(dir, file_name)`).
#' @author Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt
#' @export
#' @seealso [load_data()], [tally_repeats()], [write_specs_file()]
#' @examples
#' data_test <- data.frame(
#'   reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
#'   reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
#'   reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
#' )
#' write_files(dat = data_test, dir = tempdir())
#' run(dir = tempdir())
#'
write_files <- function(
  dat,
  dir = getwd(),
  file_dat = "data.dat",
  file_specs = "data.spc",
  minage = 0,
  maxage = NULL,
  refage = NULL,
  minusage = NULL,
  plusage = NULL,
  biasopt = NULL,
  sigopt = NULL,
  knotages = NULL
) {
  # check inputs
  if (!is.data.frame(dat)) {
    cli::cli_abort("Input 'dat' must be a data frame or tibble")
  }
  maxobs <- max(dat[, -1])
  # fill in any missing default values
  if (is.null(maxage)) {
    maxage <- ceiling(1.2 * maxobs / 5) * 5
    cli::cli_alert_info(
      "Max age not specified; using {maxage} which is the multiple of 5 which is >120% of the observed maximum"
    )
  }

  write_data_file(
    dat = dat,
    dir = dir,
    file_name = file_dat,
    minage = minage,
    maxage = maxage,
    refage = refage,
    minusage = minusage,
    plusage = plusage
  )
  write_specs_file(
    dir = dir,
    file_name = file_specs,
    nreaders = ifelse(names(dat)[1] == "count", ncol(dat) - 1, ncol(dat)),
    biasopt = biasopt,
    sigopt = sigopt,
    maxage = maxage,
    knotages = knotages
  )
}
