
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

