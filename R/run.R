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
  Output <- AgeingError::ProcessResults(
    Species = "AgeingError",
    SaveDir = directory,
    CalcEff = FALSE,
    verbose = FALSE
  )
  return(Output)
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
#' @export
#' @author Kelli F. Johnson
#' @return Invisibly return model output.
run <- function(directory, file_data, file_specs) {
  inputs <- prepare_inputs(
    file_data = file.path(directory, file_data),
    file_specs = file.path(directory, file_specs)
  )
  outputs <- prepare_run(inputs = inputs, directory = directory)
  return(invisible(outputs))
}
