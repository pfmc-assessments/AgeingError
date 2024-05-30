#' Simulate double reading data
#'
#' A function to generate simulated double reading data with given properties
#'
#' @param Nreaders The number of ageing readers
#' @param M True natural mortality
#' @param SelexForm Form of selectivity-at-age
#' (logistic selex-at-age is the only one that is implemented).
#' @param ErrorParams Error type
#' CV in the following equation:  Var[AgeRead] = (CV*TrueAge)^2
#' @param BiasParams Bias type
#' b in the following equation: E[AgeRead] = b*TrueAge
#' @param SelexParams Selectivity parameters, which
#' are standard to the logistic equation.
#' @param ReadsMat Matrix describing number of reads per reader combination.
#' Where each row specifies how many reads (in the first column)
#' have a particular pattern of double reads
#' (in the second through Nreaders+1 columns).
#' @param RecCv CV of recruitment, and it shoudl be noted that
#' recruitment is assumed to be stationary over time.
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
#' @examples
#' # Parameters for generating data
#' # This represents 2 unique readers
#' # Row 1 -- Otoliths read only once by reader
#' # Row 2 -- Otoliths read twice by reader 1
#' # Row 2 -- Otoliths read only once by reader 2
#' # Row 4 -- Otoliths read twice by reader 2
#' # Row 5 -- Otoliths read once by reader 1 and once by reader 2
#' ReadsMat <- structure(matrix(
#'   nrow = 5, ncol = 5,
#'   c(
#'     rep(25, 5),
#'     1, 1, 0, 0, 1,
#'     0, 1, 0, 0, 0,
#'     0, 0, 1, 1, 1,
#'     0, 0, 0, 1, 0
#'   )
#' ), dimnames = list(
#'   c(
#'     "Reader1_Only", "Reader1_DoubleReads",
#'     "Reader2_Only", "Reader2_DoubleReads",
#'     "Reader1_&_Reader2"
#'   ),
#'   c(
#'     "NumberOfReads",
#'     "Reader1", "Reader1_DoubleReads",
#'     "Reader2", "Reader2_DoubleReads"
#'   )
#' ))
#'
#' # Generate data
#' set.seed(2)
#' AgeReads <- SimulatorFn(
#'   Nreaders = 4, M = 0.2,
#'   SelexForm = "Logistic",
#'   SelexParams = c(5, 0.2), BiasParams = c(1, 1, 1.1, 1.1),
#'   ErrorParams = c(0.2, 0.2, 0.2, 0.2), ReadsMat = ReadsMat,
#'   RecCv = 0.6, RecAr1 = 0.8, Amax = 100
#' )
SimulatorFn <- function(
    Nreaders, M,
    SelexForm, ErrorParams, BiasParams, SelexParams,
    ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100) {
  RecDev <- stats::rnorm(Amax, mean = 0, sd = RecCv)
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
          Row <- c(
            Row,
            round(Ages[IndexI] * BiasParams[ReaderI] +
              stats::rnorm(1, mean = 0, sd = Ages[IndexI] * ErrorParams[ReaderI]))
          )
        }
      }
      AgeReads <- rbind(AgeReads, Row)
    }
  }
  return(AgeReads)
}
