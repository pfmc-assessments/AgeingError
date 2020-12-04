#' Make a row matrix
#'
#' ## The function is currently defined as
#' function (Input)
#' {
#'   if (is.vector(Input)) {
#'       Output <- t(as.matrix(Input))
#'   }
#'   if (!is.vector(Input)) {
#'       Output <- as.matrix(Input)
#'   }
#'   Output
#' }
#'
#' @param Input input to be converted into a row matrix
#' @author James T. Thorson
#' @export


rMx <-
  function(Input) {
    if (is.vector(Input)) {
      Output <- t(as.matrix(Input))
    }
    if (!is.vector(Input)) {
      Output <- as.matrix(Input)
    }
    Output
  }
