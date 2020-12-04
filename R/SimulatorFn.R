#' Simulate double reading data
#'
#' A function to generate simulated double reading data with given properties
#'
#' @param Nreaders The number of ageing readers
#' @param M True natural mortality
#' @param SelexForm Form of selectivity
#' @param ErrorParams Error type
#' @param BiasParams Bias type
#' @param SelexParams Selectivity parameters
#' @param ReadsMat Matrix describing number of reads per reader combination
#' @param RecCv CV of recruitment
#' @param RecAr1 First-order autoregressive coefficient for recruitment
#' @param Amax True maximum age
#'
#' @return Returns a simulated double read matrix
#'
#' @author James T. Thorson
#'
#' @references Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australias southern and eastern scalefish
#' and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
#'
#' @export
#'
SimulatorFn <-
  function(Nreaders, M, SelexForm, ErrorParams, BiasParams, SelexParams,
           ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100) {
    RecDev <- rnorm(Amax, mean = 0, sd = RecCv)
    for (i in 2:length(RecDev)) {
      RecDev[i] <- RecDev[i] * sqrt(1 - RecAr1) + RecDev[i - 1] * sqrt(RecAr1)
    }
    AgeStruct <- 1 * exp(-M * 1:Amax) * exp(RecDev - RecCv^2 / 2)
    if (SelexForm == "Logistic") {
      SelexAtAge <- 1 / (1 + exp((SelexParams[1] - 1:Amax) * SelexParams[2]))
    } else {
      stop("Selex not implemented")
    }
    Ages <- sample(x = 1:Amax, size = sum(ReadsMat[, 1]), prob = AgeStruct * SelexAtAge, replace = TRUE)

    AgeReads <- vector()

    IndexI <- 0
    for (RowI in 1:nrow(ReadsMat)) {
      for (OtolithI in 1:ReadsMat[RowI, 1]) {
        IndexI <- IndexI + 1
        Row <- vector()
        for (ReaderI in 1:Nreaders) {
          if (ReadsMat[RowI, ReaderI + 1] == 0) {
            Row <- c(Row, NA)
          } else {
            Row <- c(Row, round(Ages[IndexI] * BiasParams[ReaderI] + rnorm(1, mean = 0, sd = Ages[IndexI] * ErrorParams[ReaderI])))
          }
        }
        AgeReads <- rbind(AgeReads, Row)
      }
    }
    return(AgeReads)
  }
