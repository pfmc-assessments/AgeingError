#' Deprecated function renamed to load_specs()
#'
#' @param ... Any arguments associated with the deprecated function
#' @description
#' `r lifecycle::badge("deprecated")`
#' CreateSpecs() has been renamed as [load_specs()] to better clarify its purpose
#' @author Ian G. Taylor
#' @export
#' @seealso [load_specs()]
CreateSpecs <-
  function(...) {
    lifecycle::deprecate_stop(
      when = "2.1.0",
      what = "CreateSpecs()",
      with = "load_specs()"
    )
  }

#' Read the ageing error specifications
#'
#' @param SpecsFile Filename for input specifications.
#' @param DataSpecs The output from load_data()
#' @param verbose Return messages to the console (TRUE/FALSE)
#' @export
#' @author Andre E. Punt
load_specs <- function(SpecsFile = "data.spc",
                       DataSpecs,
                       verbose = FALSE) {
  MatchTable <- function(Table,
                         Char1 = NULL,
                         Char2 = NULL,
                         Char3 = NULL,
                         Char4 = NULL,
                         Char5 = NULL) {
    ii <- rep(T, length(Table[, 1]))
    if (!is.null(Char1)) {
      ii <- ii & (Table[, 1] == Char1)
    }
    if (!is.null(Char2)) {
      ii <- ii & (Table[, 2] == Char2)
    }
    if (!is.null(Char3)) {
      ii <- ii & (Table[, 3] == Char3)
    }
    if (!is.null(Char4)) {
      ii <- ii & (Table[, 4] == Char4)
    }
    if (!is.null(Char5)) {
      ii <- ii & (Table[, 5] == Char5)
    }
    ii <- seq(1:length(Table[, 1]))[ii]
    return(ii)
  }

  Specs <- read.table(
    SpecsFile,
    comment.char = "?",
    fill = TRUE,
    blank.lines.skip = TRUE,
    stringsAsFactors = FALSE,
    col.names = 1:10
  )

  # Read in the initial values
  ModelSpecs <- vector(mode = "list", length = DataSpecs$MaxReader)
  IndexA <- MatchTable(Specs, Char1 = "#", Char2 = "reader") + 1
  NumBias <- 0
  NumSigma <- 0
  for (Ireader in 1:DataSpecs$MaxReader) {
    DefaultList <- list(
      BiasOpt = NULL,
      SigOpt = NULL,
      BiasPar = NULL,
      BiasLow = NULL,
      BiasHi = NULL,
      BiasUsed = NULL,
      SigmaPar = NULL,
      SigmaLow = NULL,
      SigmaHi = NULL,
      SigmaUsed = NULL
    )
    DefaultList$BiasOpt <- as.numeric(Specs[IndexA, 2])
    DefaultList$SigOpt <- as.numeric(Specs[IndexA, 3])

    # Check for valud bias options
    if (DefaultList$BiasOpt >= 0 & !DefaultList$BiasOpt %in% c(0, 1, 2)) {
      cli::cli_alert_warning(
        "Error specifying bias option for reader {Ireader}; Bias option {DefaultList$BiasOpt} is not implemented- stopping"
      )
    }
    if (DefaultList$SigOpt >= 0 & !DefaultList$SigOpt %in% c(1:8)) {
      cli::cli_alert_warning(
        "Error specifying variance option for reader {Ireader}; Variance option {DefaultList$SigOpt} is not implemented- stopping"
      )
    }

    IndexA <- IndexA + 1
    ModelSpecs[[Ireader]] <- DefaultList
  }

  # Parameters defining Sigmas (Spline)
  Index <- MatchTable(
    Specs,
    Char1 = "#",
    Char2 = "Spline",
    Char3 = "specifications"
  )
  xvals <- matrix(0, nrow = DataSpecs$MaxReader, ncol = 100)
  nknots <- rep(0, DataSpecs$MaxReader)
  for (Ireader in 1:DataSpecs$MaxReader) {
    if (ModelSpecs[[Ireader]]$SigOpt == 5) {
      nknots[Ireader] <- as.numeric(Specs[Index + 1, 1])
      for (IDcnt in 1:nknots[Ireader]) {
        xvals[Ireader, IDcnt] <- as.numeric(Specs[Index +
          2, IDcnt])
      }
      Index <- Index + 2
      cli::cli_alert_info("Spline used to define SD for reader {Ireader}")
      cli::cli_alert_info("Selected knots are located at {paste(xvals[Ireader, 1:nknots[Ireader]], collapse = ' ')}")
    }
  }

  Index <- MatchTable(
    Specs,
    Char1 = "#",
    Char2 = "Linear",
    Char3 = "specifications"
  )
  # linear model specifications
  xvalsL <- matrix(0, nrow = DataSpecs$MaxReader, ncol = 100)
  nknotsL <- rep(0, DataSpecs$MaxReader)
  for (Ireader in 1:DataSpecs$MaxReader) {
    if (ModelSpecs[[Ireader]]$SigOpt == 6) {
      nknotsL[Ireader] <- as.numeric(Specs[Index + 1, 1])
      for (IDcnt in 1:nknotsL[Ireader]) {
        xvalsL[Ireader, IDcnt] <- as.numeric(Specs[Index +
          2, IDcnt])
      }
      Index <- Index + 2
      if (xvalsL[Ireader, 1] != DataSpecs$MinAge) {
        print("First age must be 1")
      }
      if (xvalsL[Ireader, nknotsL[Ireader]] != DataSpecs$MaxAge) {
        print("last age must be MaxAge")
      }
      cli::cli_alert_info("Linear interpolation used to define SD for reader {Ireader}")
      cli::cli_alert_info("Selected knots are located at ")
      cli::cli_alert_info("Selected knots are located at {paste(xvalsL[Ireader, 1:nknotsL[Ireader]], collapse = ' ')}")
    }
  }

  # Read in the initial values
  IndexB <- MatchTable(Specs, Char1 = "Bias_Pars") + 1
  IndexC <- MatchTable(Specs, Char1 = "Sigma_Pars") + 1
  NumBias <- 0
  NumSigma <- 0
  for (Ireader in 1:DataSpecs$MaxReader) {
    DefaultList <- list(
      BiasOpt = NULL,
      SigOpt = NULL,
      BiasPar = NULL,
      BiasLow = NULL,
      BiasHi = NULL,
      BiasUsed = NULL,
      SigmaPar = NULL,
      SigmaLow = NULL,
      SigmaHi = NULL,
      SigmaUsed = NULL
    )
    DefaultList$BiasOpt <- ModelSpecs[[Ireader]]$BiasOpt
    DefaultList$SigOpt <- ModelSpecs[[Ireader]]$SigOpt
    Nbias <- 0
    if (DefaultList$BiasOpt == 1) {
      Nbias <- 1
    }
    if (DefaultList$BiasOpt == 2) {
      Nbias <- 3
    }
    if (Nbias > 0) {
      for (Ibias in 1:Nbias) {
        DefaultList$BiasPar <- c(DefaultList$BiasPar, as.numeric(Specs[
          IndexB,
          3
        ]))
        DefaultList$BiasLow <- c(DefaultList$BiasLow, as.numeric(Specs[
          IndexB,
          1
        ]))
        DefaultList$BiasHi <- c(DefaultList$BiasHi, as.numeric(Specs[
          IndexB,
          2
        ]))
        DefaultList$BiasUsed <- c(DefaultList$BiasUsed, as.numeric(Specs[
          IndexB,
          4
        ]))
        IndexB <- IndexB + 1
      }
    }
    NumBias <- NumBias + Nbias
    Nsigma <- 0
    if (DefaultList$SigOpt == 1) {
      Nsigma <- 1
    }
    if (DefaultList$SigOpt == 2) {
      Nsigma <- 3
    }
    if (DefaultList$SigOpt == 3) {
      Nsigma <- 3
    }
    if (DefaultList$SigOpt == 5) {
      Nsigma <- nknots[Ireader]
    }
    if (DefaultList$SigOpt == 6) {
      Nsigma <- nknotsL[Ireader]
    }
    if (DefaultList$SigOpt == 7) {
      Nsigma <- 2
    }
    if (DefaultList$SigOpt == 8) {
      Nsigma <- 2
    }
    if (Nsigma > 0) {
      for (Isigma in 1:Nsigma) {
        DefaultList$SigmaPar <- c(DefaultList$SigmaPar, as.numeric(Specs[
          IndexC,
          3
        ]))
        DefaultList$SigmaLow <- c(DefaultList$SigmaLow, as.numeric(Specs[
          IndexC,
          1
        ]))
        DefaultList$SigmaHi <- c(DefaultList$SigmaHi, as.numeric(Specs[
          IndexC,
          2
        ]))
        DefaultList$SigmaUsed <- c(DefaultList$SigmaUsed, as.numeric(Specs[
          IndexC,
          4
        ]))
        IndexC <- IndexC + 1
      }
    }
    NumSigma <- NumSigma + Nsigma
    ModelSpecs[[Ireader]] <- DefaultList
  }

  Outs <- NULL
  Outs$ModelSpecs <- ModelSpecs
  Outs$NumSigma <- NumSigma
  Outs$NumBias <- NumBias
  Outs$xvalsL <- xvalsL
  Outs$nknotsL <- nknotsL
  Outs$xvals <- xvals
  Outs$nknots <- nknots
  if (verbose) {
    print(str(Outs))
  }
  return(Outs)
}
