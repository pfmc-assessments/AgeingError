#' Run ageing error model
#'
#' Run the Punt et al. (2008) ADMB-based ageing error model from within R.
#'
#' @details The premise of Punt *et al.* (2008) is to calculate the likelihood
#' of model parameters given an observed data set of otolith age reads from
#' multiple age readers. For each reader/lab, two parameters are defined, one
#' for standard deviation and one for bias. The model calculates the expected
#' age of each read and the standard deviation of a normally distributed
#' reading error given the true age of an otolith. These relationships can be
#' linear or curvilinear.
#'
#' The true age is obviously an unobserved process and can be considered a
#' random effect. Thus, the software computes the likelihood while summing
#' across all possible discrete values for the true age of each otolith. This
#' true age requires a hyperdistribution that represents the prior probability
#' that an otolith is any given age. The hyperdistribution is controlled by a
#' set of hyperparameters and the parameters that govern the standard deviation
#' and bias of each age reader/lab. Specifically, one hyperparameter is
#' estimated for every age between and including the `MinusAge` and `PlusAge`.
#' Ages outside of this range have a prior proportion at age defined as a
#' loglinear deviation from the proportion at age for the extreme ages, i.e.,
#' `MinusAge` and `PlusAge`. The slope of these loglinear deviations thus
#' constitutes an additional 1 or 2 fixed effect parameters. The true
#' proportion at age is then calculated from these fixed effects and loglinear
#' slope parameters by normalizing the resulting distribution such that it sums
#' to one.
#'
#' @param Data This is the data set with the first column being an integer
#'   providing the number of otoliths that are included in the row and the
#'   subsequent columns are the reader or lab estimated ag,e where each
#'   reader/lab has a unique reading error and bias. The modeling framework
#'   allows for, at most, 15 readers, i.e., 16 columns. There should not be any
#'   identical rows in the data frame because otoliths that have the exact same
#'   read from every reader/lab should be combined into a single row with the
#'   count as the first column. If you failed to combine identical rows prior
#'   to running the model, you will be alerted with an error and the `XXX.rep`
#'   file will have a properly formatted data which can be' cut-pasted into a
#'   `XXX.dat` file for use. Missing reads from a given reader/lab should be
#'   entered as `-999`. Order your reader/lab columns such that similar
#'   readers/labs are located next to one another because columns to the right
#'   can mirror columns to their immediate left in terms of parameter
#'   estimates.
#' @param SigOpt This a vector with one entry for each reader (i.e.,
#'   `length(SigOpt) == NCOL(Data) -1`). Each entry specifies the functional
#'   form of reading error as a function of true age. Possible entries include
#'   the following:
#'   \describe{
#'     \item{-[0-9]+}{
#'       Mirror the standard deviation of another reader, where the negative
#'       integer corresponds to the column of the reader/lab that is being
#'       mirrored minus one, e.g., `-1` causes it to mirror reader/lab 1, for
#'       which data is stored in the second column of `Data`. This number must
#'       be lower than -1 times the current position in the vector.
#'     }
#'     \item{0}{
#'       No error. But, there could be potential bias.
#'     }
#'     \item{1}{
#'       Constant coefficient of variation, i.e., a 1-parameter linear
#'       relationship of the standard deviation with the true age.
#'     }
#'     \item{2}{
#'       Curvilinear standard deviation, i.e., a 3-parameter Hollings-form
#'       relationship of standard deviation with true age.
#'     }
#'     \item{3}{
#'       Curvilinear coefficient of variation, i.e., a 3-parameter
#'       Hollings-form relationship of coefficient of variation with true age.
#'     }
#'     \item{5}{
#'       Spline with estimated slope at beginning and end where the number of
#'       parameters is 2 + number of knots.
#'     }
#'     \item{6}{
#'       Linear interpolation with a first knot of 1 and a last knot of the
#'       maximum age, i.e., `MaxAge`.
#'     }
#' }
#' @param KnotAges Ages associated with each knot. This is a necessary input
#'   for `SigOpt = 5` or `SigOpt = 6`.
#' @param BiasOpt A vector with one entry for each reader/lab specifying the
#'   type of bias specific to each reader. Positive values lead to estimated
#'   parameters and negative values are used for shared parameters between
#'   readers, just like with `SigOpt`. Parameter sharing is common when there
#'   is more than one reader in a lab working together to refine their methods
#'   such that they have matching techniques. Possible entries include the
#'   following:
#'   \describe{
#'     \item{-[0-9]+}{
#'       Mirror the bias of another reader, where the negative integer
#'       corresponds to the column of the reader/lab that is being mirrored
#'       minus one, e.g., `-1` causes it to mirror reader/lab 1, for which data
#'       is stored in the second column of `Data`. This number must be lower
#'       than -1 times the current position in the vector.
#'     }
#'     \item{0}{
#'       Unbiased, where at least one reader has to be unbiased.
#'     }
#'     \item{1}{
#'       Constant coefficient of variation, i.e., a 1-parameter linear
#'       relationship of bias with true age.
#'     }
#'     \item{2}{
#'       Curvilinear, i.e., a 2-parameter Hollings-form relationship of bias
#'       with true age.
#'     }
#'   }
#'
#'   An example entry for the situation where you have seven readers and you
#'   assume that the first reader is unbiased, readers 2-7 have a curvilinear
#'   bias, reader 3 shares parameters with reader 2, reader 5 shares parameters
#'   with reader 4, and reader 7 shares parameters with reader 6 would look
#'   like `c(0, 2, -2, 2, -4, 2, -6)`.
#' @param NDataSets This is generally `1` and other values are not implemented.
#' @param MinAge An integer, specifying the minimum possible "true" age.
#' @param MaxAge An integer, specifying the maximum possible "true" age.
#' @param RefAge An arbitrarily chosen age from which "true" age-composition
#'   fixed-effects are calculated as an offset. This has no effect on the
#'   answer but could potentially effect estimation speed.
#' @param MinusAge The minimum age for which an age-specific age-composition is
#'   estimated. Ages below `MinusAge` have "true" proportion-at-age
#'   (\eqn{P_{a}}) estimated as
#'   \deqn{P_a = P_{MinusAge}*exp^{(\beta*(MinusAge - a))}},
#'   where beta is an estimated log-linear trend in the "true"
#'   proportion-at-age. If `MinusAge` = `MinAge`, beta is not estimated.
#' @param PlusAge Identical to `MinusAge` except defining the age above with
#' age-specific age composition is not estimated.
#' @param MaxSd An upper bound on possible values for the standard deviation
#' of reading error.
#' @param MaxExpectedAge Set to MaxAge.
#' @param SaveFile Directory where `agemat.exe` is located and where all ADMB
#'   intermediate and output files should be located. If `AdmbFile` is specified
#'   then `agemat.exe` is copied from that directory to `SaveFile`.
#' @param EffSampleSize Indicating whether effective sample size should be
#'   calculated. Missing values in the data matrix will cause this to be
#'   ineffective, in which case this should be set to `0`.
#' @param Intern A logical input that controls the amount of output displayed,
#'   where `TRUE` indicates that ADMB output should be displayed in R and
#'   `FALSE` leads to the suppression of this information.
#' @param AdmbFile An optional character entry that specifies the directory
#'   from which `agemat.exe` is to be copied from to `SaveFile`.
#' @param JustWrite A logical input that allows just the data files to be
#'   written without running ADMB executable.
#' @param CallType Either `"system"` or `"shell"` depending on Operating System
#'   or how R is being run. The default is `"system"`.
#' @param ExtraArgs A string of characters providing extra arguments passed to
#'   ADMB. The default is `" -est"`.
#' @param verbose A logical input that controls the amount of feedback users
#'   receive from the program. The default is to provide the most output as
#'   possible with `verbose = TRUE`.
#' @author James T. Thorson, Ian J. Stewart, Andre E. Punt, Ian G. Taylor
#' @export
#' @seealso
#' * `StepwiseFn()` will run multiple models.
#' * `PlotOutputFn()` will help summarize the output.
#' @examples
#' example(SimulatorFn)
#' \dontrun{
#' utils::write.csv(AgeReads,
#'   file = file.path(getwd(), "Simulated_data_example.csv")
#' )
#' }
#'
#' ##### Format data
#' Nreaders <- ncol(AgeReads)
#' # Change NA to -999 (which the Punt software considers missing data)
#' AgeReads <- ifelse(is.na(AgeReads), -999, AgeReads)
#'
#' # Potentially eliminate rows that are only read once
#' # These rows have no information about reading error, but are potentially
#' # informative about latent age-structure. It is unknown whether eliminating
#' # these rows degrades estimation of error and bias, and is currently
#' # recommended to speed up computation
#' if (FALSE) {
#'   KeepRow <- ifelse(
#'     rowSums(ifelse(AgeReads == -999, 0, 1), na.rm = TRUE) <= 1,
#'     FALSE, TRUE
#'   )
#'   AgeReads <- AgeReads[KeepRow, ]
#' }
#'
#' # AgeReads2 is the correctly formatted data object
#' AgeReads2 <- rMx(c(1, AgeReads[1, ]))
#'
#' # Combine duplicate rows
#' for (RowI in 2:nrow(AgeReads)) {
#'   DupRow <- NA
#'   for (PreviousRowJ in 1:nrow(AgeReads2)) {
#'     if (all(
#'       AgeReads[RowI, 1:Nreaders] == AgeReads2[PreviousRowJ, 1:Nreaders + 1]
#'     )) {
#'       DupRow <- PreviousRowJ
#'     }
#'   }
#'   if (is.na(DupRow)) { # Add new row to AgeReads2
#'     AgeReads2 <- rbind(AgeReads2, c(1, AgeReads[RowI, ]))
#'   }
#'   if (!is.na(DupRow)) { # Increment number of samples for previous duplicate
#'     AgeReads2[DupRow, 1] <- AgeReads2[DupRow, 1] + 1
#'   }
#' }
#'
#' ######## Determine settings for ADMB
#' # Define minimum and maximum ages for integral across unobserved ages
#' MinAge <- 1
#' MaxAge <- ceiling(max(AgeReads2[, -1]) / 10) * 10
#' BiasOpt <- c(0, -1, 0, -3)
#' SigOpt <- c(1, -1, 6, -3)
#' # Necessary for SigOpt option 5 or 6
#' KnotAges <- list(NA, NA, c(1, 10, 20, MaxAge), NA)
#'
#' ##### Run the model (MAY TAKE 5-10 MINUTES)
#' \dontrun{
#' fileloc <- file.path(tempdir(), "age")
#' dir.create(fileloc, showWarnings = FALSE)
#' RunFn(
#'   Data = AgeReads2, SigOpt = SigOpt, KnotAges = KnotAges,
#'   BiasOpt = BiasOpt,
#'   NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge, RefAge = 10,
#'   MinusAge = 1, PlusAge = 30, SaveFile = fileloc,
#'   AdmbFile = file.path(system.file("executables",
#'     package = "nwfscAgeingError"
#'   ), .Platform$file.sep),
#'   EffSampleSize = 0, Intern = FALSE, JustWrite = FALSE, CallType = "shell"
#' )
#' }
RunFn <- function(Data,
                  SigOpt,
                  KnotAges,
                  BiasOpt,
                  NDataSets = 1,
                  MinAge,
                  MaxAge,
                  RefAge,
                  MinusAge,
                  PlusAge,
                  MaxSd,
                  MaxExpectedAge,
                  SaveFile,
                  EffSampleSize = 0,
                  Intern = TRUE,
                  AdmbFile = NULL,
                  JustWrite = FALSE,
                  CallType = "system",
                  ExtraArgs = " -est",
                  verbose = TRUE) {
  # add slash to end of directories so that nobody has to waste as much time
  # debugging as Ian just did
  SaveFile <- paste0(SaveFile, "/")

  # Copy ADMB file
  if (!is.null(AdmbFile)) {
    AdmbFile <- paste0(AdmbFile, "/")
    # Check for missing file before trying to copy
    if (is.na(file.info(file.path(AdmbFile, "agemat.exe"))$size)) {
      warning(
        "executable 'agemat.exe' not found in\n",
        AdmbFile
      )
    }
    if (verbose) {
      cat(
        "copying 'agemat.exe' from\n", AdmbFile,
        "\nto\n", SaveFile, "\n"
      )
    }
    file.copy(
      from = file.path(AdmbFile, "agemat.exe"),
      to = file.path(SaveFile, "agemat.exe"), overwrite = TRUE
    )
  }
  # Check for missing file
  if (is.na(file.info(file.path(SaveFile, "agemat.exe"))$size)) {
    stop(
      "executable 'agemat.exe' not found in\n",
      SaveFile
    )
  }

  # Check for errors
  Nreaders <- ncol(Data) - 1
  for (ReaderI in 1:Nreaders) {
    if ((SigOpt[ReaderI] == 5 | SigOpt[ReaderI] == 6) &
      is.na(KnotAges[[ReaderI]][1])) {
      stop("Must specify KnotAges for any reader with SigOpt 5 or 6")
    }
  }

  # Check for specification errors
  for (ReaderI in 1:Nreaders) {
    if ((SigOpt[ReaderI] < 0 & SigOpt[ReaderI] <= (-ReaderI)) |
      (BiasOpt[ReaderI] < 0 & BiasOpt[ReaderI] <= (-ReaderI))) {
      stop("Mirrored readers must mirror a lower numbered reader")
    }
  }

  # Write DAT file
  datfile <- file.path(SaveFile, "agemat.dat")
  # simple functions to avoid repeated code
  writeLine <- function(x, file = datfile, append = TRUE) {
    write(x = x, file = file, append = append)
  }
  writeTable <- function(x, file = datfile, append = TRUE) {
    utils::write.table(
      x = x, file = file, append = append,
      row.names = FALSE, col.names = FALSE
    )
  }

  # first line creates new file (not appended to existing file)
  writeLine(c("# Maximum number of readers", Nreaders), append = FALSE)
  # subsequent lines use default append = TRUE
  writeLine(c("# Number of data sets", NDataSets))
  writeLine(c("# Number of points per data set", nrow(Data)))
  writeLine(c("# Readers per data set", ncol(Data) - 1))
  writeLine("# Which readers per data set")
  writeTable(rMx(1:(ncol(Data) - 1)))
  writeLine(c("# Minimum age", MinAge))
  writeLine(c("# Maximum age", MaxAge))
  writeLine(c("# Reference age", RefAge))
  writeLine(c("# Minus groups", MinusAge))
  writeLine(c("# Plus groups", PlusAge))
  # write table of options for bias
  writeLine("# Option for bias")
  writeTable(rMx(BiasOpt))
  # write table of options for SD
  writeLine("# Option for standard deviation")
  writeTable(rMx(SigOpt))
  writeLine(c("# Option for effective sample size", EffSampleSize))
  writeLine(c("# Use Par File (1=Yes)", 0))
  # Write knots related to splines
  writeLine("\n# Number and location of knots for any splines")
  for (ReaderI in 1:Nreaders) {
    if (SigOpt[ReaderI] == 5 | SigOpt[ReaderI] == 6) {
      writeTable(rMx(c(length(KnotAges[[ReaderI]]), KnotAges[[ReaderI]])))
    }
  }
  # Write initial values
  # Bias
  writeLine("\n# Min, Max, Init, Phase for Bias")
  for (BiasI in 1:Nreaders) {
    # No bias
    if (BiasOpt[BiasI] <= 0) {}
    # Linear bias
    if (BiasOpt[BiasI] == 1) {
      writeTable(rMx(c(0.001, 3, 1, 2)))
    }
    # Curvilinear bias = 0.5+Par1 + (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)))*(1.0-mfexp(-Par2*(float(Age1)-1)))
    # Starting value must be non-zero
    if (BiasOpt[BiasI] == 2) {
      writeTable(rMx(c(0.001, 10, 1, 2)))
      writeTable(rMx(c(-10, 1, 0.01, 2)))
      writeTable(rMx(c(0.001, MaxAge * 2, MaxAge, 2)))
    }
  }
  # Sigma
  writeLine("\n# Min, Max, Init, Phase for Sigma")
  for (SigI in 1:Nreaders) {
    # No error
    if (SigOpt[SigI] <= 0) {}
    # Linear CV
    if (SigOpt[SigI] == 1) {
      writeTable(rMx(c(0.001, 3, 0.1, 2)))
    }
    # Curvilinear SD = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if (SigOpt[SigI] == 2) {
      writeTable(rMx(c(0.001, 100, 1, 2)))
      writeTable(rMx(c(-10, 1, 0.01, 2)))
      writeTable(rMx(c(0.001, 100, 10, 2)))
    }
    # Curvilinear CV
    # Curvilinear CV = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if (SigOpt[SigI] == 3) {
      writeTable(rMx(c(0.001, 3, 0.1, 2)))
      writeTable(rMx(c(-10, 1, 0.01, 2)))
      writeTable(rMx(c(0.001, 3, 0.1, 2)))
    }
    # Spline with estimated derivative at beginning and end
    # (Params 1-N: knot parameters; N+1 and N+2: derivative at beginning and end)
    if (SigOpt[SigI] == 5) {
      for (ParI in 1:(2 + length(KnotAges[[SigI]]))) {
        writeTable(rMx(c(-10.0, 10.0, 1.0, 1)))
      }
    }
    # Spline with derivative at beginning and end fixed at zero
    # (Params 1-N: knot parameters)
    if (SigOpt[SigI] == 6) {
      for (ParI in 1:length(KnotAges[[SigI]])) {
        writeTable(rMx(c(-10.0, 10.0, 1.0, 1)))
      }
    }
  }
  # Probs (i.e. age-composition probability relative to reference age)
  writeLine("\n# Min, Max, Phase for Probs")
  writeTable(rMx(c(-20, 20, 2)))
  # Slopes
  writeLine("\n# Min, Max, Init, Phase for slopes")
  for (DataSetI in 1:NDataSets) {
    if (MaxAge > PlusAge) {
      writeTable(rMx(c(-10, 0, 0, 1)))
    }
    if (MinAge < MinusAge) {
      writeTable(rMx(c(-10, 0, 0, 1)))
    }
  }
  # Write dataset
  writeLine("\n# Data set")
  writeTable(Data)
  writeLine(c("# Test number", 123456))

  # Run ADMB file
  if (JustWrite == FALSE) {
    setwd(SaveFile)
    if (CallType == "shell") {
      # This may need to have the location pasted onto it depending upon file structure
      Output <- shell(paste0("agemat.exe", ExtraArgs), intern = Intern)
    }
    if (CallType == "system") {
      Output <- system(paste0("agemat.exe", ExtraArgs), intern = Intern)
    }
    # Admb = scan(paste(SaveFile,"agemat.par",sep=""),comment.char="#",quiet=TRUE)
  }
  # return(Output)
}
