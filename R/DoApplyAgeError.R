#' Run the ageing error optimization routine
#'
#' @param Species A string that will be used to create file names. Typically,
#'   users will use the common name for the species of interest, especially if
#'   you are saving files from multiple species in a single directory. Though,
#'   the default is `"AgeingError"`.
#' @param DataSpecs A data object returned from [load_data()].
#' @param ModelSpecsInp A specification object returned from [load_specs()].
#' @param AprobWght,SlopeWght Numeric values passed to the model. The defaults
#'   are 1e-06 and 0.01. Andre originally had these hard coded from his
#'   workspace. TODO: decide if they should be passed in the specifications or
#'   data files.
#' @param SaveDir A path, relative or absolute, to a directory where the
#'   results will be saved. The directory need not exist currently as it will
#'   be created dynamically.
#' @param verbose A logical specifying if messages should be printed. The
#'   default is to **NOT** print, i.e., `verbose = FALSE`.
#'
#' @export
#' @author Andre E. Punt
DoApplyAgeError <- function(Species = "AgeingError",
                            DataSpecs,
                            ModelSpecsInp,
                            AprobWght = 1e-06,
                            SlopeWght = 0.01,
                            SaveDir = getwd(),
                            verbose = FALSE) {
  fs::dir_create(SaveDir)
  SaveFile <- file.path(
    SaveDir,
    paste0(Species, ".lda")
  )
  if (verbose) {
    cli::cli_inform(
      "Results will be saved to {SaveFile}"
    )
  }

  # Extract material from the data specs
  MinAge <- DataSpecs$MinAge
  MaxAge <- DataSpecs$MaxAge
  MaxReader <- DataSpecs$MaxReader
  NDataSet <- DataSpecs$NDataSet
  Npnt <- DataSpecs$Npnt
  ReadPnt <- DataSpecs$ReadPnt
  TheData <- DataSpecs$TheData
  Nread <- DataSpecs$NReaders
  MinusA <- DataSpecs$MinusA
  PlusA <- DataSpecs$PlusA
  RefAge <- DataSpecs$RefAge

  # extract the Sigma and bias options
  BiasOpt <- rep(NA, MaxReader)
  SigOpt <- rep(NA, MaxReader)
  NumBias <- ModelSpecsInp$NumBias
  NumSig <- ModelSpecsInp$NumSigma
  Bias_LO <- NULL
  Bias_HI <- NULL
  Bias_INIT <- NULL
  Sigma_LO <- NULL
  Sigma_HI <- NULL
  Sigma_INIT <- NULL
  ModelSpecs <- ModelSpecsInp$ModelSpecs
  if (NumBias > 0) {
    BiasParMap <- as.factor(NULL)
  } else {
    BiasParMap <- factor(NA)
  }
  if (NumSig > 0) {
    SDParMap <- as.factor(NULL)
  } else {
    SDParMap <- factor(NA)
  }
  IBiasCnt <- 0
  ISigCnt <- 0
  for (Ireader in 1:MaxReader) {
    TheSpecs <- ModelSpecs[[Ireader]]
    SigOpt[Ireader] <- TheSpecs$SigOpt
    BiasOpt[Ireader] <- TheSpecs$BiasOpt
    if (BiasOpt[Ireader] > 0) {
      for (II in 1:length(TheSpecs$BiasLow)) {
        IBiasCnt <- IBiasCnt + 1
        if (TheSpecs$BiasUsed[II] > 0) {
          Bias_LO <- c(Bias_LO, TheSpecs$BiasLow[II])
          Bias_HI <- c(Bias_HI, TheSpecs$BiasHi[II])
          BiasParMap <- c(BiasParMap, factor(IBiasCnt))
        } else {
          BiasParMap <- c(BiasParMap, factor(NA))
        }
      }
      Bias_INIT <- c(Bias_INIT, TheSpecs$BiasPar)
    }
    if (SigOpt[Ireader] > 0) {
      for (II in 1:length(TheSpecs$SigmaLow)) {
        ISigCnt <- ISigCnt + 1
        if (TheSpecs$SigmaUsed[II] > 0) {
          Sigma_LO <- c(Sigma_LO, TheSpecs$SigmaLow[II])
          Sigma_HI <- c(Sigma_HI, TheSpecs$SigmaHi[II])
          SDParMap <- c(SDParMap, factor(ISigCnt))
        } else {
          SDParMap <- c(SDParMap, factor(NA))
        }
      }
      Sigma_INIT <- c(Sigma_INIT, TheSpecs$SigmaPar)
    }
  }

  # Find number of reads by data set and if any readers are perfect
  # Number of readers for this dataset
  Nreads <- rep(0, NDataSet)
  # Number of data sets for which the answer is not known
  NDataSetWithoutPerfect <- 0
  # one if there is no perfect reading in this dataset
  Iperfect <- rep(0, NDataSet)
  for (IDataS in 1:NDataSet) {
    for (II in 1:Nread[IDataS]) {
      if (ReadPnt[IDataS, II] > 0) {
        Nreads[IDataS] <- Nreads[IDataS] + 1
        if (SigOpt[ReadPnt[IDataS, II]] == 4) {
          Iperfect[IDataS] <- II
        }
      }
    }
    if (Iperfect[IDataS] == 0) {
      NDataSetWithoutPerfect <- NDataSetWithoutPerfect + 1
    }
  }

  # Determine the number of parameters that should be estimated
  Nprobs <- 0
  Nslops <- 0
  for (IDataS in 1:NDataSet) {
    if (Iperfect[IDataS] == 0) {
      Nprobs <- Nprobs + (PlusA[IDataS] - MinusA[IDataS])
      if (MaxAge > PlusA[IDataS]) {
        Nslops <- Nslops + 1
      }
      if (MinAge < MinusA[IDataS]) {
        Nslops <- Nslops + 1
      }
    }
  }

  # Initial values
  Slope_LO <- rep(-10, Nslops)
  Slope_HI <- rep(1, Nslops)
  Slope_INIT <- rep(0, Nslops)

  # initial values
  Prob_LO <- rep(-20, Nprobs)
  Prob_HI <- rep(20, Nprobs)
  Prob_INIT <- rep(0, Nprobs)
  Probs <- rep(0, Nprobs)


  # Initialize the proportion parameters
  Jpnt <- 0
  for (IDataS in 1:NDataSet) {
    if (Iperfect[IDataS] == 0) {
      AgFreq <- rep(0, MaxAge + 1)
      for (Ipnt in 1:Npnt[IDataS]) {
        for (Iread in 1:Nread[IDataS]) {
          if (TheData[IDataS, Ipnt, Iread + 1] >= 0) {
            AgFreq[TheData[IDataS, Ipnt, Iread + 1] + 1] <- AgFreq[TheData[
              IDataS,
              Ipnt, Iread + 1
            ] + 1] + TheData[IDataS, Ipnt, 1]
          }
        }
      }
      Normal <- AgFreq[RefAge[IDataS] + 1]
      for (Age in (MinusA[IDataS]):(PlusA[IDataS])) {
        AgFreq[Age + 1] <- log((AgFreq[Age + 1] + 0.1) / Normal)
      }

      for (Age in MinusA[IDataS]:(RefAge[IDataS] - 1)) {
        Jpnt <- Jpnt + 1
        Probs[Jpnt] <- AgFreq[Age + 1]
      }
      for (Age in (RefAge[IDataS] + 1):(PlusA[IDataS])) {
        Jpnt <- Jpnt + 1
        Probs[Jpnt] <- AgFreq[Age + 1]
      }
    }
  }
  Prob_INIT <- Probs

  # Now apply
  data <- list(
    NDataSet = NDataSet,
    MinAge = MinAge,
    MaxAge = MaxAge,
    BiasOpt = BiasOpt,
    SigOpt = SigOpt,
    MaxReader = MaxReader,
    MinusA = MinusA,
    PlusA = PlusA,
    RefAge = RefAge,
    Iperfect = Iperfect,
    TheData = TheData,
    Npnt = Npnt,
    Nread = Nread,
    MaxNpnt = max(Npnt),
    ReadPnt = ReadPnt,
    ReaderSumm = DataSpecs$ReaderSumm,
    ReaderStruc = DataSpecs$ReaderStruc,
    MaxCells = DataSpecs$MaxCells,
    TotalN = DataSpecs$TotalN,
    EffN = DataSpecs$EffN,
    xvals = ModelSpecsInp$xvals,
    nknots = ModelSpecsInp$nknots,
    xvalsL = ModelSpecsInp$xvalsL,
    nknotsL = ModelSpecsInp$nknotsL,
    AprobWght = AprobWght,
    SlopeWght = SlopeWght
  )
  if (verbose) {
    print(str(data))
  }

  Bias_INIT_Use <- Bias_INIT
  if (is.null(Bias_INIT)) {
    Bias_INIT_Use <- 1
  }
  parameters <- list(
    Dummy = 0,
    BiasPar = Bias_INIT_Use,
    SDPar = Sigma_INIT,
    Slope = Slope_INIT,
    Probs = Prob_INIT
  )
  if (verbose) {
    print(str(parameters))
  }
  map <- list(
    BiasPar = rep(factor(NA), NumBias), SDPar = rep(factor(NA), NumSig),
    Slope = rep(factor(NA), Nslops), Probs = rep(factor(NA), Nprobs)
  )
  map <- list(Dummy = factor(NA), BiasPar = BiasParMap, SDPar = SDParMap)

  if (is.null(Bias_INIT)) {
    map <- list(Dummy = factor(NA), BiasPar = BiasParMap, SDPar = SDParMap)
  }
  if (verbose) {
    print(map)
  }

  ################## Sigma_LO <- c(0,0.01,0,0,0) Sigma_HI <- c(1,1,2,1,2)
  upper <- c(Bias_HI, Sigma_HI, Slope_HI, Prob_HI)
  lower <- c(Bias_LO, Sigma_LO, Slope_LO, Prob_LO)
  model <- TMB::MakeADFun(
    data,
    parameters,
    map = map,
    silent = TRUE,
    DLL = "AgeingError"
  )
  model$fn_orig <- model$fn

  # Debugging track
  model$fn <- function(x) {
    vv <- model$fn_orig(x)
    print(vv)
    return(vv)
  }

  # Find the model and iterate until convergence
  model <- minimizer(model, method = "both", lower, upper)
  best <- 1e+20
  print(model$fitv, digits = 10)
  while (abs(best - model$fitv) > 1e-10) {
    print("looping")
    print(best)
    best <- model$fitv
    model <- minimizer(model, method = "both", lower, upper)
    print("Objective fn:")
    print(model$fitv, digits = 10)
    cat("Difference", best - model$fitv, "\n")
  }
  model <- minimizer(model, method = "both", lower, upper, verbose = verbose)
  print(model$gr(model$env$last.par.best))
  print(model$env$last.par.best)

  # Save tesults
  SaveAll <- NULL
  SaveAll$model <- model
  SaveAll$data <- data
  SaveAll$parameters <- parameters
  SaveAll$report <- model$report()
  save(SaveAll, file = SaveFile)
  rep <- TMB::sdreport(model)
  if (verbose) {
    print(summary(rep))
  }
  SaveAll$gradient <- rep$gradient.fixed
  SaveAll$sdreport <- rep
  save(SaveAll, file = SaveFile)
  return(model)
}
