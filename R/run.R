prepare_inputs <- function(file_data, file_specs) {
  the_data <- CreateData(
    DataFile = file_data,
    NDataSet = determine_n_sets(file_data)
  )
  the_specs <- CreateSpecs(
    SpecsFile = file_specs,
    DataSpecs = the_data
  )
  return(list(
    data = the_data,
    specs = the_specs
  ))
}

prepare_run <- function(inputs, directory) {
  model <- AgeingError:::DoApplyAgeError(
    Species = "AgeingError",
    DataSpecs = inputs[["data"]],
    ModelSpecs = inputs[["specs"]],
    AprobWght = 0.000001,
    SlopeWght = 0.01,
    SaveDir = directory,
    verbose = FALSE
  )
  Output <- AgeingError:::ProcessResults(
    Species = "AgeingError",
    SaveDir = directory,
    CalcEff = FALSE,
    verbose = FALSE
  )
  return(model)
}

#' Run ageing error routine
#'
#' A wrapper for running a TMB model to estimate ageing error for a given data
#' set and specification file.
#'
#' @param file_data A string specifying the file path to a data file.
#' @param file_specs A string specifying the file path to the specifications
#'   file.
#' @param directory A string specifying a file path to a directory where you
#'   would like to save the results.
#'
#' @export
#' @author Kelli F. Johnson
#' @return Invisibly return model output.
run <- function(file_data, file_specs, directory) {
  inputs <- prepare_inputs(file_data = file_data, file_specs = file_specs)
  outputs <- prepare_run(inputs = inputs, directory = directory)
  return(invisible(outputs))
}
