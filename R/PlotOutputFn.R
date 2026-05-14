#' Deprecated function replaced by plot_output()
#'
#' @param ... Any arguments associated with the deprecated function
#' @description
#' `r lifecycle::badge("deprecated")`
#' PlotOutputFn() has been replaced by [plot_output()]
#' @author James T. Thorson, Ian G. Taylor
#' @export
#' @seealso [plot_output()]
PlotOutputFn <-
  function(...) {
    lifecycle::deprecate_stop(
      when = "2.1.1",
      what = "PlotOutputFn()",
      with = "plot_output()"
    )
  }