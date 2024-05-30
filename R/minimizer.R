#' Minimize the negative log likelihood
#'
#' Minimize the negative log likelihood using `"nlmimb"` and/or `"optim"`.
#'
#' @param model A model to be optimized.
#' @param method A string specifying the desired method to be used for the
#'   optimization routine. The options are listed in the function call, where
#'   the default is to use `"optim"`. Using both routines is an option, via
#'   `"both"`, and will lead to first optimizing the model using `"nlminb"`
#'   and then re-optimization of the model with `"optim"`. Note that when using
#'   [stats::optim()], the `"L-BFGS-B"` method is used rather than the default
#'   method of `"Nelder-Mead"`.
#' @param lower,upper Vectors of parameter bounds of the same length as the
#'   number of parameters in the model.
#' @inheritParams DoApplyAgeError
#' @export
#' @author Andre E. Punt
minimizer <- function(model,
                      method = c("optim", "nlmimb", "both"),
                      lower,
                      upper,
                      verbose = FALSE) {
  method <- match.arg(method)
  # Check parameters 'work'
  if (length(lower) > 0 & length(model$par) != length(lower)) {
    cli::cli_abort("wrong number of lower bounds")
  }
  if (length(upper) > 0 & length(model$par) != length(upper)) {
    cli::cli_abort("wrong number of upper bounds")
  }

  # Find the model and iterate until convergence
  if (method == "both" || method == "nlmimb") {
    fit <- nlminb(
      model$par,
      model$fn,
      model$gr,
      upper = upper,
      lower = lower,
      control = list(
        eval.max = 10000,
        iter.max = 10000,
        rel.tol = 1e-15,
        abs.tol = 1e-15
      )
    )
    model$par <- model$env$last.par.best
    model$fitv <- fit$objective
    if (verbose) {
      cli::cli_inform(model$fitv)
    }
  }
  if (method == "both" || method == "optim") {
    fit <- optim(
      model$par,
      model$fn,
      method = "L-BFGS-B",
      model$gr,
      upper = upper,
      lower = lower,
      control = list(maxit = 10000, factr = 1e-15)
    )
    model$par <- model$env$last.par.best
    model$fitv <- fit$value
    if (verbose) {
      cli::cli_inform(model$fitv)
    }
  }
  if (verbose) {
    print(fit)
  }

  return(model)
}
