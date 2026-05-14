#' Deprecated function replaced by run()
#'
#' @param ... Any arguments associated with the deprecated function
#' @description
#' `r lifecycle::badge("deprecated")`
#' RunFn() has been replaced by [run()]
#' @author James T. Thorson, Ian J. Stewart, Andre E. Punt, Ian G. Taylor
#' @export
#' @seealso [run()]
RunFn <-
  function(...) {
    lifecycle::deprecate_stop(
      when = "2.1.1",
      what = "RunFn()",
      with = "run()"
    )
  }