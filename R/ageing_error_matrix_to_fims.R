#' Format an ageing-error probability matrix for FIMS
#'
#' Convert a matrix with true ages in rows and observed ages in columns to the
#' long format expected for `ageing_error` input by FIMS.
#'
#' @param probabilities Numeric square matrix of probabilities. Rows represent
#'   true ages and columns represent observed ages.
#' @param ages Age values corresponding to the rows and columns of
#'   `probabilities`.
#' @param fleet Fleet name, or `NA` to apply to every fleet with age-composition
#'   data.
#' @param timing Model year, or `NA` to use as the default for every year.
#' @return A data frame with columns required for FIMS `ageing_error` input.
#' @export
ageing_error_matrix_to_fims <- function(
  probabilities,
  ages,
  fleet = NA_character_,
  timing = NA_real_
) {
  if (!is.matrix(probabilities) || !is.numeric(probabilities)) {
    cli::cli_abort("{.arg probabilities} must be a numeric matrix.")
  }
  if (nrow(probabilities) != ncol(probabilities) ||
    length(ages) != nrow(probabilities) || length(ages) == 0) {
    cli::cli_abort(
      "{.arg probabilities} must be square and have one value in {.arg ages} per row."
    )
  }
  if (!is.numeric(ages) || anyNA(ages) || any(!is.finite(ages)) ||
    any(ages != floor(ages)) || anyDuplicated(ages)) {
    cli::cli_abort("{.arg ages} must contain unique finite whole numbers.")
  }
  if (length(fleet) != 1 || (!is.character(fleet) && !is.na(fleet))) {
    cli::cli_abort("{.arg fleet} must be one character value or {.code NA}.")
  }
  if (length(timing) != 1 || !is.numeric(timing) ||
    (!is.na(timing) && (!is.finite(timing) || timing != floor(timing)))) {
    cli::cli_abort("{.arg timing} must be one whole-number year or {.code NA}.")
  }
  if (anyNA(probabilities) || any(!is.finite(probabilities)) ||
    any(probabilities < 0 | probabilities > 1)) {
    cli::cli_abort("{.arg probabilities} must contain finite values between 0 and 1.")
  }
  if (any(abs(rowSums(probabilities) - 1) > 1e-3)) {
    cli::cli_abort("Each true-age row in {.arg probabilities} must sum to 1.")
  }

  if (is.na(fleet)) {
    fleet <- NA_character_
  }

  n_ages <- length(ages)
  data.frame(
    type = "ageing_error",
    fleet = fleet,
    age = rep(ages, times = n_ages),
    length = NA_real_,
    timing = timing,
    observed = as.vector(t(probabilities)),
    unit = "proportion",
    uncertainty = as.character(rep(ages, each = n_ages))
  )
}