#' Run ageing error model
#'
#' Run the Punt et al. (2008) ADMB-based ageing error model from within R
#'
#' @param Data This is the data set as previously formatted. If the data has
#' multiple rows with identical reads, this will cause an error and the
#' "XXX.rep" file will have a properly formatted data matrix which can be
#' cut-pasted into a "XXX.dat" file for use.
#' @param SigOpt This a vector with one entry for each reader (i.e. Ncol-1
#' entries). Each entry specifies the functional form of reading error as
#' a function of true age. Possible entries include:
#' \itemize{
#'   \item{"-1", "-2", "-3", etc: This will make this reader mirror the
#'         estimated SD from another reader to it's left. "-1" causes it to
#'         mirror the estimated SD for the first reader, etc. This number has
#'         to be lower than the current entry number.}
#'   \item{"0": No error (but potentially bias)}
#'   \item{"1": Constant CV, i.e., a 1 parameter linear relationship of SD with
#'         true age.}
#'   \item{"2": Curvilinear SD, i.e., a 3 parameter Hollings-form relationship
#'         of SD with true age}
#'   \item{"3": Curvilinear with CV, i.e., a 3-parameter Hollings-form
#'         relationship of CV with true age}
#'   \item{"5": Spline with estimated slope at beginning and end (Number
#'         of params = 2 + number of knots)}
#'   \item{"6": Linear interpolation (1st knot must be 1 and last knot must
#'         be MaxAge)}
#' }
#' @param KnotAges Ages associated with (necessary for options 5 or 6)
#' @param BiasOpt A vector with one entry for each reader specifying the type
#' of bias specific to each reader. Positive values lead to estimated
#' parameters and negative values are used for shared parameters
#' between readers. Parameter sharing is common when there is more
#' than one reader in a lab working together to refine their methods
#' such that they have matching techniques. The following types of bias
#' options are available:
#' \describe{
#'   \item{0}{Unbiased}
#'   \item{1}{Constant CV: a 1-parameter linear relationship of bias
#'     with true age}
#'   \item{2}{Curvilinear: a 2-parameter Hollings-form relationship
#'     of bias with true age}
#' }
#' An example entry for the situation where you have seven readers and you assume
#' that the first reader is unbiased, readers 2-7 have a curvilinear
#' bias, reader 3 shares parameters with reader 2, reader 5 shares parameters
#' with reader 4, and reader 7 shares parameters with reader 6 would look like
#' \code{c(0, 2, -2, 2, -4, 2, -6)}.
#' @param NDataSets This is generally "1" and other values are not implemented
#' in the current R-code.
#' @param MinAge The minimum possible "true" age
#' @param MaxAge The maximum possible "true" age
#' @param RefAge An arbitrarily chosen age from which "true" age-composition
#' fixed-effects are calculated as an offset. This has no effect on the answer,
#' but could potentially effect estimation speed.
#' @param MinusAge The minimum age for which an age-specific age-composition is
#' estimated. Ages below this MinusAge have "true" proportion-at-age (P_a)
#' estimated as P_a = P_MinusAge*exp(beta*(MinusAge - a)), where beta is an
#' estimated log-linear trend in the "true" proportion-at-age.
#' If MinusAge = MinAge, beta is not estimated.
#' @param PlusAge Identical to MinusAge except defining the age above with
#' age-specific age-composition is not estimated.
#' @param MaxSd An upper bound on possible values for the standard deviation
#' of reading error
#' @param MaxExpectedAge Set to MaxAge
#' @param SaveFile Directory where "agemat.exe" is located and where all ADMB
#' intermediate and output files should be located. If AdmbFile is specified
#' then "agemat.exe" is copied from that directory to SaveFile
#' @param EffSampleSize Indicating whether effective sample size should be
#' calculated. Missing values in the data matrix will cause this to be
#' ineffective, in which case this should be set to "0"
#' @param Intern "TRUE" indicates that ADMB output should be displayed in R;
#' "FALSE" does not.
#' @param AdmbFile Optional directory from which "agemat.exe" is to be copied
#' to SaveFile
#' @param JustWrite Switch to allow data files to be written without running
#' ADMB executable.
#' @param CallType Either "system" or "shell" depending on Operating System
#' or how R is being run.
#' @param ExtraArgs Extra arguments passed to ADMB. Default is " -est".
#' @param verbose Provide more feedback about function progress?
#' @author James T. Thorson, Ian J. Stewart, Andre E. Punt, Ian G. Taylor
#' @export
#' @seealso \code{\link{StepwiseFn}}, \code{\link{PlotOutputFn}}

RunFn <-
  function(Data, SigOpt, KnotAges, BiasOpt, NDataSets, MinAge, MaxAge, RefAge,
           MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile,
           EffSampleSize = 0, Intern = TRUE, AdmbFile = NULL, JustWrite = FALSE,
           CallType = "system", ExtraArgs = " -est", verbose = TRUE)
{

  # add slash to end of directories so that nobody has to waste as much time
  # debugging as Ian just did
  SaveFile <- paste0(SaveFile, "/")

  # Copy ADMB file
  if(!is.null(AdmbFile)){
    AdmbFile <- paste0(AdmbFile, "/")
    # Check for missing file before trying to copy
    if(is.na(file.info(file.path(AdmbFile,"agemat.exe"))$size)){
      warning("executable 'agemat.exe' not found in\n",
              AdmbFile)
    }
    if(verbose){
      cat("copying 'agemat.exe' from\n", AdmbFile,
          "\nto\n", SaveFile,'\n')
    }
    file.copy(from = file.path(AdmbFile,"agemat.exe"),
              to = file.path(SaveFile,"agemat.exe"), overwrite = TRUE)
  }
  # Check for missing file
  if(is.na(file.info(file.path(SaveFile,"agemat.exe"))$size)){
    stop("executable 'agemat.exe' not found in\n",
         SaveFile)
  }

  # Check for errors
  Nreaders = ncol(Data)-1
  for(ReaderI in 1:Nreaders){
    if( (SigOpt[ReaderI] == 5 | SigOpt[ReaderI] == 6) &
       is.na(KnotAges[[ReaderI]][1]) ) {
      stop("Must specify KnotAges for any reader with SigOpt 5 or 6")
    }
  }

  # Check for specification errors
  for(ReaderI in 1:Nreaders){
    if( (SigOpt[ReaderI]<0 & SigOpt[ReaderI]<=(-ReaderI)) |
       (BiasOpt[ReaderI]<0 & BiasOpt[ReaderI]<=(-ReaderI)) ){
      stop("Mirrored readers must mirror a lower numbered reader")
    }
  }

  # Write DAT file
  datfile <- file.path(SaveFile, "agemat.dat")
  # simple functions to avoid repeated code
  writeLine <- function(x, file = datfile, append=TRUE){
    write(x = x, file = file, append = append)
  }
  writeTable <- function(x, file = datfile, append=TRUE){
    write.table(x = x, file = file, append = append,
                row.names = FALSE, col.names = FALSE)
  }

  # first line creates new file (not appended to existing file)
  writeLine(c("# Maximum number of readers",Nreaders), append = FALSE)
  # subsequent lines use default append = TRUE
  writeLine(c("# Number of data sets",NDataSets))
  writeLine(c("# Number of points per data set",nrow(Data)))
  writeLine(c("# Readers per data set",ncol(Data)-1))
  writeLine("# Which readers per data set")
  writeTable(rMx(1:(ncol(Data)-1)))
  writeLine(c("# Minimum age",MinAge))
  writeLine(c("# Maximum age",MaxAge))
  writeLine(c("# Reference age",RefAge))
  writeLine(c("# Minus groups",MinusAge))
  writeLine(c("# Plus groups",PlusAge))
  # write table of options for bias
  writeLine("# Option for bias")
  writeTable(rMx(BiasOpt))
  # write table of options for SD
  writeLine("# Option for standard deviation")
  writeTable(rMx(SigOpt))
  writeLine(c("# Option for effective sample size",EffSampleSize))
  writeLine(c("# Use Par File (1=Yes)",0))
  # Write knots related to splines
  writeLine("\n# Number and location of knots for any splines")
  for(ReaderI in 1:Nreaders){
    if(SigOpt[ReaderI] == 5 | SigOpt[ReaderI] == 6){
      writeTable( rMx(c(length(KnotAges[[ReaderI]]),KnotAges[[ReaderI]])))
    }
  }
  # Write initial values
  # Bias
  writeLine("\n# Min, Max, Init, Phase for Bias")
  for(BiasI in 1:Nreaders){
    # No bias
    if(BiasOpt[BiasI]<=0){}
    # Linear bias
    if(BiasOpt[BiasI] == 1){
      writeTable(rMx(c(0.001,3,1,2)))
    }
    # Curvilinear bias = 0.5+Par1 + (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)))*(1.0-mfexp(-Par2*(float(Age1)-1)))
    # Starting value must be non-zero
    if(BiasOpt[BiasI] == 2){
      writeTable(rMx(c(0.001,10,1,2)))
      writeTable(rMx(c(-10,1,0.01,2)))
      writeTable(rMx(c(0.001,MaxAge*2,MaxAge,2)))
    }
  }
  # Sigma
  writeLine("\n# Min, Max, Init, Phase for Sigma")
  for(SigI in 1:Nreaders){
    # No error
    if(SigOpt[SigI]<=0){}
    # Linear CV
    if(SigOpt[SigI] == 1){
      writeTable(rMx(c(0.001,3,0.1,2)))
    }
    # Curvilinear SD = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if(SigOpt[SigI] == 2){
      writeTable(rMx(c(0.001,100,1,2)))
      writeTable(rMx(c(-10,1,0.01,2)))
      writeTable(rMx(c(0.001,100,10,2)))
    }
    # Curvilinear CV
    # Curvilinear CV = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
    # Starting value must be non-zero
    if(SigOpt[SigI] == 3){
      writeTable(rMx(c(0.001,3,0.1,2)))
      writeTable(rMx(c(-10,1,0.01,2)))
      writeTable(rMx(c(0.001,3,0.1,2)))
    }
    # Spline with estimated derivative at beginning and end
    # (Params 1-N: knot parameters; N+1 and N+2: derivative at beginning and end)
    if(SigOpt[SigI] == 5){
      for(ParI in 1:(2+length(KnotAges[[SigI]]))){
        writeTable(rMx(c(-10.0,10.0,1.0,1)))
      }
    }
    # Spline with derivative at beginning and end fixed at zero
    # (Params 1-N: knot parameters)
    if(SigOpt[SigI] == 6){
      for(ParI in 1:length(KnotAges[[SigI]])){
        writeTable(rMx(c(-10.0,10.0,1.0,1)))
      }
    }
  }
  # Probs (i.e. age-composition probability relative to reference age)
  writeLine("\n# Min, Max, Phase for Probs")
  writeTable(rMx(c(-20,20,2)))
  # Slopes
  writeLine("\n# Min, Max, Init, Phase for slopes")
  for(DataSetI in 1:NDataSets){
    if(MaxAge > PlusAge){
      writeTable(rMx(c(-10,0,0,1)))
    }
    if(MinAge < MinusAge){
      writeTable(rMx(c(-10,0,0,1)))
    }
  }
  # Write dataset
  writeLine("\n# Data set")
  writeTable(Data)
  writeLine(c("# Test number",123456))

  # Run ADMB file
  if(JustWrite == FALSE){
    setwd(SaveFile)
    if(CallType == "shell"){
      # This may need to have the location pasted onto it depending upon file structure
      Output = shell( paste0("agemat.exe",ExtraArgs),intern=Intern)
    }
    if(CallType == "system"){
      Output = system( paste0("agemat.exe",ExtraArgs),intern=Intern)
    }
    #Admb = scan(paste(SaveFile,"agemat.par",sep=""),comment.char="#",quiet=TRUE)
  }
  #return(Output)
}
