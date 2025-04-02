prepare_inputs <- function(file_data, file_specs) {
  the_data <- load_data(
    DataFile = file_data,
    NDataSet = determine_n_sets(file_data)
  )
  the_specs <- load_specs(
    SpecsFile = file_specs,
    DataSpecs = the_data
  )
  return(list(
    data = the_data,
    specs = the_specs
  ))
}

prepare_run <- function(inputs, directory) {
  model <- AgeingError::DoApplyAgeError(
    Species = "AgeingError",
    DataSpecs = inputs[["data"]],
    ModelSpecs = inputs[["specs"]],
    AprobWght = 0.000001,
    SlopeWght = 0.01,
    SaveDir = directory,
    verbose = FALSE
  )
  output <- AgeingError::ProcessResults(
    Species = "AgeingError",
    SaveDir = directory,
    CalcEff = FALSE,
    verbose = FALSE
  )
  return(list(model = model, output = output))
}

#' Run ageing error routine
#'
#' A wrapper for running a TMB model to estimate ageing error for a given data
#' set and specification file.
#'
#' @param directory A string specifying a file path to a directory where you
#'   would like to save the results.
#' @param file_data A string specifying the data file within 'directory'.
#' @param file_specs A string specifying the specifications file within 'directory'.
#'
#' @seealso [write_files()]
#' @examples
#' \dontrun{
#' data_test <- data.frame(
#'   reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
#'   reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
#'   reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
#' )
#' write_files(dat = data_test, dir = tempdir())
#' out <- run(dir = tempdir())
#' # see estimated parameters
#' out$model$par
#' # see model selection results
#' out$output$ModelSelection
#' # see ageing error matrices
#' out$output$ErrorAndBiasArray
#' # add to an SS3 model (assumes the model already has a single
#' # ageing error matrix and a maxage <= maxage in the ageing error model)
#' ss3_inputs <- r4ss::SS_read()
#' maxage <- ss3_inputs$dat$Nages
#' ss3_inputs$dat$ageerror <- out$output$ErrorAndBiasArray[c("Expected_age", "SD"), 1 + 0:maxage, "Reader 1"] |>
#'   as.data.frame()
#' r4ss::SS_write(inputlist = ss3_inputs)
#' }
#' @export
#' @author Kelli F. Johnson
#' @return Invisibly return model output.
run <- function(directory, file_data = "data.dat", file_specs = "data.spc") {
  # Check if the directory exists
  if (!dir.exists(directory)) {
    cli::cli_abort("The specified directory does not exist.")
  }
  # Check if the data file exists
  if (!file.exists(file.path(directory, file_data))) {
    cli::cli_abort("The specified data file does not exist: {file.path(directory, file_data)}")
  }
  # Check if the specs file exists
  if (!file.exists(file.path(directory, file_specs))) {
    cli::cli_abort("The specified specifications file does not exist: {file.path(directory, file_specs)}")
  }
  inputs <- prepare_inputs(
    file_data = file.path(directory, file_data),
    file_specs = file.path(directory, file_specs)
  )
  outputs <- prepare_run(inputs = inputs, directory = directory)
  return(invisible(outputs))
}
