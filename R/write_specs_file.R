#' Write a specifications file for the AgeingError package
#'
#' The resulting file has the format required by the [load_specs()] function.
#' In the future, this step could be bypassed by creating the list output
#' from [load_specs()] directly.
#' This function is based on a subset of the RunFn() function in the older ADMB
#' version of this package.
#'
#' @param dir Directory where the specifications file will be saved.
#' @param file_name Name of the specifications file.
#' @param nreaders An integer, specifying the number of readers.
#' @param biasopt A vector with one entry for each reader specifying the
#'   type of bias specific to each reader. Positive values lead to estimated
#'   parameters and negative values are used for shared parameters between
#'   readers. Parameter sharing (mirroring) is common when there
#'   is more than one reader in a lab working together to refine their methods
#'   such that they have matching techniques.
#'   If NULL is passed, the default is `rep(0, nreaders)` which specifies that
#'   all readers are unbiased.
#'
#'   Possible entries include the following:
#'   \describe{
#'     \item{\code{-[0-9]+}}{
#'       Mirror the bias of another reader, where the negative integer
#'       corresponds to the column of the reader that is being mirrored
#'       minus one, e.g., `-1` causes it to mirror reader 1. Only lower-numbered
#'       readers can be mirrored.
#'     }
#'     \item{0}{
#'       Unbiased, where at least one reader has to be unbiased.
#'     }
#'     \item{1}{
#'       Constant coefficient of variation, i.e., a 1-parameter linear
#'       relationship of bias with true age.
#'     }
#'     \item{2}{
#'       Curvilinear, i.e., a 2-parameter Hollings-form relationship of bias
#'       with true age.
#'     }
#'   }
#'
#'   An example entry for the situation where you have seven readers and you
#'   assume that the first reader is unbiased, readers 2-7 have a curvilinear
#'   bias, reader 3 shares parameters with reader 2, reader 5 shares parameters
#'   with reader 4, and reader 7 shares parameters with reader 6 would look
#'   like `c(0, 2, -2, 2, -4, 2, -6)`.
#' @param sigopt A vector with one entry for each reader.
#'   Each entry specifies the functional
#'   form of reading error as a function of true age. Positive values lead to
#'   estimated parameters and negative values are used for shared parameters
#'   between readers.
#'   If NULL is passed, the default is `c(1, rep(-1, nreaders - 1))` which
#'   specifies a constant CV in ageing error which is shared among all readers.
#'
#'   Possible entries include the following:
#'   \describe{
#'     \item{\code{-[0-9]+}}{
#'       Mirror the standard deviation of another reader, where the negative
#'       integer corresponds to the column of the reader that is being
#'       mirrored minus one, e.g., `-1` causes it to mirror reader 1, for
#'       which data is stored in the second column of `Data`. This number must
#'       be lower than -1 times the current position in the vector.
#'     }
#'     \item{0}{
#'       No error. But, there could be potential bias.
#'     }
#'     \item{1}{
#'       Constant coefficient of variation, i.e., a 1-parameter linear
#'       relationship of the standard deviation with the true age.
#'     }
#'     \item{2}{
#'       Curvilinear standard deviation, i.e., a 3-parameter Hollings-form
#'       relationship of standard deviation with true age.
#'     }
#'     \item{3}{
#'       Curvilinear coefficient of variation, i.e., a 3-parameter
#'       Hollings-form relationship of coefficient of variation with true age.
#'     }
#'     \item{5}{
#'       Spline with estimated slope at beginning and end where the number of
#'       parameters is 2 + number of knots. Supported when `knotages` is provided.
#'     }
#'     \item{6}{
#'       Linear interpolation with a first knot of 1 and a last knot of the
#'       maximum age, i.e., `MaxAge`. Supported when `knotages` is provided.
#'     }
#'     \item{7}{
#'       A linear change in the standard deviation of random age-reading error,
#'       $\sigma_a$, with age. This option has two parameters that need to be
#'       specified for each pair of independent readers in the specifications
#'       file.
#'     }
#'     \item{8}{
#'        A linear change in the coefficient of variation of random age-reading
#'        error, $CV_a$ , with age. This option has two parameters that need to
#'        be specified for each pair of independent readers in the
#'        specifications file.
#'     }
#'   }
#' @param knotages A list of knot ages for each reader. This is required when
#'   `sigopt = 5` or `sigopt = 6` and must have one element per reader.
#' @param maxage The maximum possible "true" age.
#'   This is only required in the specifications file if `biasopt = 2`.
#' @return Invisibly returns the path to the specifications file (`file.path(dir, file_name)`).
#' @author Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt
#' @export
#' @seealso [write_files()], [write_specs_file()], [load_specs()]
#' @examples
#' data_test <- data.frame(
#'   reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
#'   reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
#'   reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
#' )
#' data_file <- write_data_file(data_test, dir = tempdir(), file_name = "test.dat")
#' specs_file <- write_specs_file(dir = tempdir(), nreaders = 3, file_name = "test.spc")
#' data <- load_data(DataFile = data_file)
#' specs <- load_specs(SpecsFile = specs_file, DataSpecs = data)
#'
write_specs_file <- function(
    dir = getwd(),
    file_name = "data.spc",
    nreaders,
    biasopt = NULL,
    sigopt = NULL,
    knotages = NULL,
    maxage) {
  # check inputs
  if (!dir.exists(dir)) {
    cli::cli_alert_info("Directory does not exist; creating it")
    dir.create(dir)
  }
  if (!is.character(file_name)) {
    cli::cli_abort("Input 'file_name' must be a character string")
  }
  if (!is.numeric(nreaders)) {
    cli::cli_abort("Input 'nreaders' must be a numeric value")
  }
  if (!is.null(biasopt) && (!is.numeric(biasopt) || length(biasopt) != nreaders)) {
    cli::cli_abort("Input 'biasopt' must be a numeric vector of length 'nreaders'")
  }
  if (!is.null(sigopt) && (!is.numeric(sigopt) || length(sigopt) != nreaders)) {
    cli::cli_abort("Input 'sigopt' must be a numeric vector of length 'nreaders'")
  }

  # fill in any missing default values
  if (is.null(biasopt)) {
    biasopt <- rep(0, nreaders)
    cli::cli_alert_info(
      "'biasopt' not specified; settings all readers to unbiased"
    )
  }
  if (is.null(sigopt)) {
    sigopt <- c(1, rep(-1, nreaders - 1))
    cli::cli_alert_info(
      "'sigopt' not specified; settings all readers to share a constant CV parameter"
    )
  }
  if (any(sigopt %in% 5:6)) {
    if (is.null(knotages)) {
      cli::cli_abort("'knotages' must be specified when 'sigopt' includes 5 or 6")
    }
    if (!is.list(knotages) || length(knotages) != nreaders) {
      cli::cli_abort("'knotages' must be a list with one element per reader")
    }
    for (ireader in which(sigopt %in% 5:6)) {
      if (!is.numeric(knotages[[ireader]]) || length(knotages[[ireader]]) < 2) {
        cli::cli_abort(
          "knotages[[{ireader}]] must be a numeric vector with at least two knot ages for sigopt {sigopt[ireader]}"
        )
      }
    }
  }

  # create reader specs table
  reader_specs <- data.frame(reader = 1:nreaders, biasopt, sigopt)

  # bias parameters
  biaspars <- NULL
  # fill in parameter lines (defaults below based on old RunFn() function)
  for (ireader in 1:nreaders) {
    # No bias
    if (biasopt[ireader] <= 0) {
      newpars <- NULL
    }
    # Linear bias
    if (biasopt[ireader] == 1) {
      newpars <- data.frame(low = 0.001, high = 3, init = 1, on_off = 1)
    }
    # Curvilinear bias = 0.5+Par1 + (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)))*(1.0-mfexp(-Par2*(float(Age1)-1)))
    # Starting value must be non-zero
    if (biasopt[ireader] == 2) {
      newpars <- data.frame(
        low = c(0.001, -10, 0.001),
        high = c(10, 1, maxage * 2),
        init = c(1, 0.01, maxage),
        on_off = 1
      )
    }
    biaspars <- rbind(biaspars, newpars)
  }

  # sigma parameters
  sigpars <- NULL
  for (ireader in 1:nreaders) {
    # No error
    if (sigopt[ireader] <= 0) {
      newpars <- NULL
    }
    # Linear CV
    if (sigopt[ireader] == 1) {
      newpars <- data.frame(low = 0.001, high = 3, init = 0.1, on_off = 1)
    }
    # Curvilinear SD = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if (sigopt[ireader] == 2) {
      newpars <- data.frame(
        low = c(0.001, -10, 0.001),
        high = c(100, 1, 100),
        init = c(1, 0.01, 10),
        on_off = 1
      )
    }
    # Curvilinear CV
    # Curvilinear CV = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if (sigopt[ireader] == 3) {
      newpars <- data.frame(
        low = c(0.001, -10, 0.001),
        high = c(3, 1, 3),
        init = c(0.1, 0.01, 0.1),
        on_off = 1
      )
    }
    # Spline with estimated derivative at beginning and end
    if (sigopt[ireader] == 5) {
      newpars <- data.frame(
        low = rep(-10, length(knotages[[ireader]])),
        high = rep(40, length(knotages[[ireader]])),
        init = rep(0, length(knotages[[ireader]])),
        on_off = 1
      )
    }
    # Piecewise linear with specified knots
    if (sigopt[ireader] == 6) {
      newpars <- data.frame(
        low = rep(-10, length(knotages[[ireader]])),
        high = rep(40, length(knotages[[ireader]])),
        init = rep(0, length(knotages[[ireader]])),
        on_off = 1
      )
    }
    # Linear change in SD with age
    if (sigopt[ireader] == 7) {
      newpars <- data.frame(
        low = c(0, 0),
        high = c(1, 2),
        init = c(0.2, 0.5),
        on_off = 1
      )
    }
    # Linear change in CV with age
    if (sigopt[ireader] == 8) {
      newpars <- data.frame(
        low = c(0, 0),
        high = c(1, 2),
        init = c(0.2, 0.5),
        on_off = 1
      )
    }
    sigpars <- rbind(sigpars, newpars)
  }

  # write file
  cli::cli_alert_info(
    "Writing specifications file to {file.path(dir, file_name)}"
  )
  # write reader specifications
  names(reader_specs)[1] <- "# reader" # require format for load_specs()
  readr::write_delim(
    reader_specs,
    file = file.path(dir, file_name),
    quote = "none"
  )

  # write knot specifications for spline/linear sigma options
  if (any(sigopt == 5)) {
    readr::write_lines(
      c(" ", "# Spline specifications"),
      file = file.path(dir, file_name),
      append = TRUE
    )
    for (ireader in which(sigopt == 5)) {
      readr::write_lines(
        as.character(length(knotages[[ireader]])),
        file = file.path(dir, file_name),
        append = TRUE
      )
      readr::write_lines(
        paste(knotages[[ireader]], collapse = " "),
        file = file.path(dir, file_name),
        append = TRUE
      )
    }
    readr::write_lines(
      " ",
      file = file.path(dir, file_name),
      append = TRUE
    )
  }
  if (any(sigopt == 6)) {
    readr::write_lines(
      c(" ", "# Linear specifications"),
      file = file.path(dir, file_name),
      append = TRUE
    )
    for (ireader in which(sigopt == 6)) {
      readr::write_lines(
        as.character(length(knotages[[ireader]])),
        file = file.path(dir, file_name),
        append = TRUE
      )
      readr::write_lines(
        paste(knotages[[ireader]], collapse = " "),
        file = file.path(dir, file_name),
        append = TRUE
      )
    }
    readr::write_lines(
      " ",
      file = file.path(dir, file_name),
      append = TRUE
    )
  }

  # write bias parameters
  readr::write_lines(
    c(" ", "Bias_Pars (low high init, on/off)"),
    file = file.path(dir, file_name),
    append = TRUE
  )
  if (!is.null(biaspars)) {
    names(biaspars)[1] <- paste("#", names(biaspars)[1])
    readr::write_delim(
      biaspars,
      col_names = FALSE,
      file = file.path(dir, file_name),
      quote = "none",
      append = TRUE
    )
  }

  # write sigma parameters
  readr::write_lines(
    c(" ", "Sigma_Pars (low high init, on/off)"),
    file = file.path(dir, file_name),
    append = TRUE
  )
  if (!is.null(sigpars)) {
    names(sigpars)[1] <- paste("#", names(sigpars)[1])
    readr::write_delim(
      sigpars,
      col_names = FALSE,
      file = file.path(dir, file_name),
      quote = "none",
      append = TRUE
    )
  }

  return(invisible(file.path(dir, file_name)))
}
