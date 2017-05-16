#' Stepwise model selection
#'
#' Code for running stepwise model selection
#'
#' @param SearchMat               Description needed
#' @param Data                    Description needed
#' @param NDataSets               Description needed
#' @param MinAge                  Description needed
#' @param MaxAge                  Description needed
#' @param RefAge                  Description needed
#' @param MaxSd                   Description needed
#' @param MaxExpectedAge          Description needed
#' @param SaveFile                Description needed
#' @param EffSampleSize           Description needed
#' @param Intern                  Description needed
#' @param InformationCriterion    Description needed
#'
#' @references Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australias southern and eastern scalefish
#' and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
#'
#' @author James T. Thorson
#'
#' @export
#'
#' @examples
#'   \dontrun{
#'  # Run initial model
#'  example(nwfscAgeingError)
#'
#'  ####################
#'  #
#'  # Stepwise selection
#'  #
#'  ####################
#'
#'  # Parameters
#'  MaxAge = ceiling(max(AgeReads2)/10)*10
#'  MinAge = 1
#'
#'  ##### Stepwise selection
#'  StartMinusAge = 1
#'  StartPlusAge = 30
#'
#'  # Define matrix explaining stepwise model selection options
#'  # One row for each reader + 2 rows for PlusAge (the age where the proportion-at-age begins to decrease exponentially with increasing age) and MinusAge (the age where the proportion-at-age begins to decrease exponentially with decreasing age)
#'  # Each element of a given row is a possible value to search across for that reader
#'  SearchMat = array(NA, dim=c(Nreaders*2+2,7), dimnames=list(c(paste("Error_Reader",1:Nreaders),paste("Bias_Reader",1:Nreaders),"MinusAge","PlusAge"),paste("Option",1:7)))
#'    # Readers 1 and 3 search across options 1-3 for ERROR
#'    SearchMat[c(1,3),1:3] = rep(1,2) %o% c(1,2,3)
#'    # Reader 2 mirrors reader 1
#'    SearchMat[2,1] = -1
#'    # Reader 4 mirrors reader 3
#'    SearchMat[4,1] = -3
#'    # Reader 1 has no BIAS
#'    SearchMat[5,1] = 0
#'    # Reader 2 mirrors reader 1
#'    SearchMat[6,1] = -1
#'    # Reader 3 search across options 0-2 for BIAS
#'    SearchMat[7,1:3] = c(1,2,0)
#'    # Reader 4 mirrors reader 3
#'    SearchMat[8,1] = -3
#'    # MinusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
#'    SearchMat[9,1:7] = c(StartMinusAge,StartMinusAge-10,StartMinusAge-4,StartMinusAge-1,StartMinusAge+1,StartMinusAge+4,StartMinusAge+10)
#'      SearchMat[9,1:7] = ifelse(SearchMat[9,1:7]<MinAge,NA,SearchMat[9,1:7])
#'    # PlusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
#'    SearchMat[10,1:7] = c(StartPlusAge,StartPlusAge-10,StartPlusAge-4,StartPlusAge-1,StartPlusAge+1,StartPlusAge+4,StartPlusAge+10)
#'      SearchMat[10,1:7] = ifelse(SearchMat[10,1:7]>MaxAge,NA,SearchMat[10,1:7])
#'
#'  # Run model selection
#'  # This outputs a series of files
#'    # 1. "Stepwise - Model loop X.txt" -- Shows the AIC/BIC/AICc value for all different combinations of parameters arising from changing one parameter at a time according to SearchMat during loop X
#'    # 2. "Stepwise - Record.txt" -- The Xth row of IcRecord shows the record of the Information Criterion for all trials in loop X, while the Xth row of StateRecord shows the current selected values for all parameters at the end of loop X
#'    # 3. Standard plots for each loop
#'  # WARNING: One run of this stepwise model building example can take 8+ hours, and should be run overnight
#'  StepwiseFn(SearchMat=SearchMat, Data=AgeReads2, NDataSets=1, MinAge=MinAge, MaxAge=MaxAge, RefAge=10, MaxSd=40, MaxExpectedAge=MaxAge+10, SaveFile=DateFile, InformationCriterion=c("AIC","AICc","BIC")[3])
#' }

StepwiseFn <-
function(SearchMat, Data, NDataSets, MinAge, MaxAge, RefAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize=0, Intern=TRUE, InformationCriterion="AIC", SelectAges=TRUE){

  # Define variables
  Nages = MaxAge+1
  Nreaders = ncol(Data)-1
  ParamVecOpt = SearchMat[,1]
  Stop = FALSE
  IcRecord = NULL
  StateRecord = NULL
  OuterIndex = 0

  # Continue searching until Stop==TRUE
  while(Stop==FALSE){

    # Increment and intialize variables
    OuterIndex = OuterIndex+1
    Index = 0
    IcVec = NULL
    ParamMat = NULL
    Rep = NULL
    ParamVecOptPreviouslyEstimates = FALSE

    # Loop across all combinations of parameters obtained from one change in the current parameters
    for(VarI in 1:nrow(SearchMat)){
    for(ValueI in 1:length(na.omit(SearchMat[VarI,]))){

      # Update the current vector of parameters
      ParamVecCurrent = ParamVecOpt
      ParamVecCurrent[VarI] = na.omit(SearchMat[VarI,])[ValueI]

      # Decide if this ParamVecCurrent should be run
      if( all(ParamVecCurrent==ParamVecOpt)&ParamVecOptPreviouslyEstimates==FALSE | !all(ParamVecCurrent==ParamVecOpt)){

        # If running the current optimum, change so that it won't run again this loop
        if(all(ParamVecCurrent==ParamVecOpt)) ParamVecOptPreviouslyEstimates=TRUE

        # Make a new file for ADMB
        RunFile = paste(SaveFile,"Run\\",sep="")
        dir.create(RunFile)
        file.copy(from=paste(SaveFile,"agemat.exe",sep=""), to=paste(RunFile,"agemat.exe",sep=""))

        # Increment Index
        Index = Index+1
        print(paste("Loop=",OuterIndex," Run=",Index," StartTime=",date(),sep=""))

        # Configure and run model
        SigOpt = ParamVecCurrent[1:Nreaders]
        BiasOpt = ParamVecCurrent[Nreaders+1:Nreaders]
        MinusAge = ParamVecCurrent[2*Nreaders+1]
        PlusAge = ParamVecCurrent[2*Nreaders+2]
        RunFn(SigOpt=SigOpt, KnotAges=KnotAges, BiasOpt=BiasOpt, Data=Data, NDataSets=NDataSets, MinAge=MinAge, MaxAge=MaxAge, RefAge=RefAge, MinusAge=MinusAge, PlusAge=PlusAge, MaxSd=MaxSd, MaxExpectedAge=MaxExpectedAge, SaveFile=RunFile)

        # Compute information criteria
        Df = as.numeric(scan(paste(RunFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[6])
        Nll = as.numeric(scan(paste(RunFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[11])
        n = sum(ifelse(Data[,-1]==-999,0,1))
        Aic = 2*Nll + 2*Df
        Aicc = Aic + 2*Df*(Df+1)/(n-Df-1)
        Bic = 2*Nll + Df*log(n)
        if(InformationCriterion=="AIC") IcVec = c(IcVec,  Aic)
        if(InformationCriterion=="AICc") IcVec = c(IcVec,  Aicc)
        if(InformationCriterion=="BIC") IcVec = c(IcVec,  Bic)
        ParamMat = rbind(ParamMat, ParamVecCurrent)
        write.table(cbind(IcVec,ParamMat),paste(SaveFile,"Stepwise - Model loop ",OuterIndex,".txt",sep=""),sep="\t",row.names=FALSE)

        # Input misclassification matrices
        Rep[[Index]] = readLines(paste(RunFile,"agemat.rep",sep=""))

     } # End if-statement for only running ParamVecOpt once per loop

    }}  # End loop accross VarI and ValueI

    IcRecord = rbind(IcRecord, IcVec)
    StateRecord = rbind(StateRecord, ParamVecOpt)
    capture.output(list(IcRecord=IcRecord, StateRecord=StateRecord), file=paste(SaveFile,"Stepwise - Record.txt",sep=""))

    # Change current vector of optimum parameters
    Max = which.max(IcVec)
    if(all(ParamMat[Max,]==ParamVecOpt)) Stop=TRUE
    ParamVecOpt = ParamMat[Max,]

    # Change boundaries for MinusAge parameter
    CurrentMinusAge = ParamVecOpt[length(ParamVecOpt)-1]
    if(SelectAges==TRUE){
      SearchMat[length(ParamVecOpt)-1,1:7] = c(CurrentMinusAge,CurrentMinusAge-10,CurrentMinusAge-4,CurrentMinusAge-1,CurrentMinusAge+1,CurrentMinusAge+4,CurrentMinusAge+10)
      SearchMat[length(ParamVecOpt)-1,1:7] = ifelse(SearchMat[length(ParamVecOpt)-1,1:7]<MinAge,NA,SearchMat[length(ParamVecOpt)-1,1:7])
    }else{
      SearchMat[length(ParamVecOpt)-1,1:7] = c( CurrentMinusAge, rep(NA,6) )
    }

    # Change boundaries for PlusAge parameter
    CurrentPlusAge = ParamVecOpt[length(ParamVecOpt)]
    if(SelectAges==TRUE){
      SearchMat[length(ParamVecOpt),1:7] = c(CurrentPlusAge,CurrentPlusAge-10,CurrentPlusAge-4,CurrentPlusAge-1,CurrentPlusAge+1,CurrentPlusAge+4,CurrentPlusAge+10)
      SearchMat[length(ParamVecOpt),1:7] = ifelse(SearchMat[length(ParamVecOpt),1:7]>MaxAge,NA,SearchMat[length(ParamVecOpt),1:7])
    }else{
      SearchMat[length(ParamVecOpt),1:7] = c( CurrentPlusAge, rep(NA,6) )
    }
    # Save image for each while loop
    writeLines(Rep[[Max]], con=paste(SaveFile,"agemat.rep",sep=""))
    PlotOutputFn(Data=Data, MaxAge=MaxAge, SaveFile=SaveFile, PlotType="JPG")

  } # End while statement

} # End StepwiseFn

