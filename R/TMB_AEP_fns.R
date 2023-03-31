Minimzer <- function(model,method="optim",lower,upper, verbose = FALSE)
{
  # Check parameters "work" 
  if (length(lower) > 0 & length(model$par) != length(lower)) { print("wrong number of lower bounds"); AA }
  if (length(upper) > 0 & length(model$par) != length(upper)) { print("wrong number of upper bounds"); AA }
  
  # Find the model and iterate until convergence
  if (method=="both" || method =="nlmimb")
   {  
    fit   <- nlminb(model$par, model$fn, model$gr,upper=upper,lower=lower,control=list(eval.max=10000,iter.max=10000,rel.tol=1.0e-15,abs.tol=1.0e-15))
    model$par <- model$env$last.par.best
    model$fitv <- fit$objective
    print(model$fitv,digits=10)
   }
  if (method=="both" || method =="optim")
   {  
    fit   <- optim(model$par, model$fn, method="L-BFGS-B",model$gr,upper=upper,lower=lower,control=list(maxit=10000,factr=1.0e-15))
    model$par <- model$env$last.par.best
    model$fitv <- fit$value
    print(model$fitv,digits=10)
  }
 if (verbose) print(fit)
  
 return(model)  
}


DoApplyAgeError <- function(Species,DataSpecs,ModelSpecsInp,SaveDir="Final", verbose = FALSE)
 {

  if (!dir.exists(SaveDir)) dir.create(SaveDir)
  
  # Save
  SaveFile <- paste(SaveDir,"/",Species,".lda",sep="")
  print(SaveFile)

  # Extract material from the data specis
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
  BiasOpt <- rep(NA,MaxReader)
  SigOpt  <- rep(NA,MaxReader)
  NumBias <- ModelSpecsInp$NumBias; NumSig <- ModelSpecsInp$NumSigma
  Bias_LO <- NULL; Bias_HI <- NULL; Bias_INIT <- NULL
  Sigma_LO <- NULL; Sigma_HI <- NULL; Sigma_INIT <- NULL
  ModelSpecs <- ModelSpecsInp$ModelSpecs
  if (NumBias>0) BiasParMap=as.factor(NULL) else BiasParMap <- factor(NA)
  if (NumSig>0) SDParMap=as.factor(NULL) else SDParMap <- factor(NA)
  IBiasCnt <- 0; ISigCnt <- 0
  for (Ireader in 1:MaxReader)
   {
    TheSpecs <- ModelSpecs[[Ireader]]
    SigOpt[Ireader] <- TheSpecs$SigOpt
    BiasOpt[Ireader] <- TheSpecs$BiasOpt
    if (BiasOpt[Ireader]>0) 
     {
      for (II in 1:length(TheSpecs$BiasLow))
       {
        IBiasCnt <- IBiasCnt+ 1 
        if (TheSpecs$BiasUsed[II]>0) 
         {
          Bias_LO <- c(Bias_LO,TheSpecs$BiasLow[II])
          Bias_HI <- c(Bias_HI,TheSpecs$BiasHi[II])
          BiasParMap <- c(BiasParMap,factor(IBiasCnt))
        }
        else
         BiasParMap <- c(BiasParMap,factor(NA))
      }  
      Bias_INIT <- c(Bias_INIT,TheSpecs$BiasPar);
     }
    if (SigOpt[Ireader]>0) 
     {
      for (II in 1:length(TheSpecs$SigmaLow))
       {
        ISigCnt <- ISigCnt+ 1 
         if (TheSpecs$SigmaUsed[II]>0) 
         {
          Sigma_LO <- c(Sigma_LO,TheSpecs$SigmaLow[II])
          Sigma_HI <- c(Sigma_HI,TheSpecs$SigmaHi[II])
          SDParMap <- c(SDParMap,factor(ISigCnt))
        }
        else
         SDParMap <- c(SDParMap,factor(NA))
      }
      Sigma_INIT <- c(Sigma_INIT,TheSpecs$SigmaPar);
     }      
   }  

  # Find number of reads by data set and if any readers are perfect
  Nreads <- rep(0,NDataSet)                               # Number of readers for this dataset
  NDataSetWithoutPerfect <-0                              # Number of data sets for which the answer is not known
  Iperfect <- rep(0,NDataSet)                             # one if there is no perfect reading in this dataset
  for (IDataS in 1:NDataSet)
   { 
    for ( II in 1:Nread[IDataS])
     if (ReadPnt[IDataS,II] > 0)
       { Nreads[IDataS] <- Nreads[IDataS] + 1;if (SigOpt[ReadPnt[IDataS,II]] == 4) Iperfect[IDataS] <- II; }
    if (Iperfect[IDataS] == 0) NDataSetWithoutPerfect <- NDataSetWithoutPerfect  + 1
   }

  # Determine the number of parameters that should be estimated
  Nprobs <- 0; Nslops <- 0; 
  for (IDataS in 1:NDataSet)
   if (Iperfect[IDataS] == 0)
    {
     Nprobs <- Nprobs + (PlusA[IDataS] - MinusA[IDataS]);
     if (MaxAge > PlusA[IDataS]) Nslops <- Nslops + 1;
     if (MinAge < MinusA[IDataS]) Nslops <- Nslops + 1;
   }
  
  # Initial values
  Slope_LO <- rep(-10,Nslops)
  Slope_HI <- rep(1,Nslops)
  Slope_INIT <- rep(0,Nslops)
  
  # initial values
  Prob_LO <- rep(-20,Nprobs)
  Prob_HI <- rep( 20,Nprobs)
  Prob_INIT <- rep(0,Nprobs)
  Probs <- rep(0,Nprobs)
  
  
  # Initialize the proportion parameters
  Jpnt <- 0;
  for (IDataS in 1:NDataSet)
   if (Iperfect[IDataS] == 0)
    {
     AgFreq <- rep(0,MaxAge+1);
     for (Ipnt in 1:Npnt[IDataS])
      for (Iread in 1:Nread[IDataS])
       if (TheData[IDataS,Ipnt,Iread+1] >= 0)
         AgFreq[TheData[IDataS,Ipnt,Iread+1]+1] <- AgFreq[TheData[IDataS,Ipnt,Iread+1]+1] + TheData[IDataS,Ipnt,1];
     Normal <- AgFreq[RefAge[IDataS]+1];
     for (Age in (MinusA[IDataS]):(PlusA[IDataS])) AgFreq[Age+1] = log( (AgFreq[Age+1]+0.1)/Normal)

     for (Age in MinusA[IDataS]:(RefAge[IDataS]-1)) { Jpnt <- Jpnt+1; Probs[Jpnt] = AgFreq[Age+1]; }
     for (Age in (RefAge[IDataS]+1):(PlusA[IDataS]))  { Jpnt <- Jpnt+1; Probs[Jpnt] = AgFreq[Age+1]; }
    }
  Prob_INIT <- Probs

  # ==========================================================================================================
  # Now apply
  # ==========================================================================================================

  data <- list(NDataSet=NDataSet,MinAge=MinAge,MaxAge=MaxAge,BiasOpt=BiasOpt,SigOpt=SigOpt,MaxReader=MaxReader,
               MinusA=MinusA,PlusA=PlusA,RefAge=RefAge,Iperfect=Iperfect,TheData=TheData,Npnt=Npnt,Nread=Nread,
               MaxNpnt=max(Npnt),ReadPnt=ReadPnt,ReaderSumm=DataSpecs$ReaderSumm,ReaderStruc=DataSpecs$ReaderStruc,
               MaxCells=DataSpecs$MaxCells,TotalN=DataSpecs$TotalN,EffN=DataSpecs$EffN,
               xvals=ModelSpecsInp$xvals,nknots=ModelSpecsInp$nknots,xvalsL=ModelSpecsInp$xvalsL,nknotsL=ModelSpecsInp$nknotsL,
               AprobWght=AprobWght,SlopeWght=SlopeWght)
  if (verbose) print(str(data))

  Bias_INIT_Use <- Bias_INIT
  if (is.null(Bias_INIT)) Bias_INIT_Use <- 1
  parameters=list(Dummy=0,BiasPar=Bias_INIT_Use,SDPar=Sigma_INIT,Slope=Slope_INIT,Probs=Prob_INIT)
  if (verbose) print(str(parameters))
  map <- list(BiasPar=rep(factor(NA),NumBias),SDPar=rep(factor(NA),NumSig),Slope=rep(factor(NA),Nslops),Probs=rep(factor(NA),Nprobs))
  map <- list(Dummy=factor(NA),BiasPar=BiasParMap,SDPar=SDParMap)
  
  if (is.null(Bias_INIT)) map <- list(Dummy=factor(NA),BiasPar=BiasParMap,SDPar=SDParMap)
  if (verbose) print(map)

  ##################
  #Sigma_LO <- c(0,0.01,0,0,0)
  #Sigma_HI <- c(1,1,2,1,2)
  upper <- c(Bias_HI,Sigma_HI,Slope_HI,Prob_HI)
  lower <- c(Bias_LO,Sigma_LO,Slope_LO,Prob_LO)
  model <- TMB::MakeADFun(data,parameters,map=map,silent=T,DLL='AgeingError')
  model$fn_orig <- model$fn

  # Debugging track
  model$fn <- function(x) { vv <- model$fn_orig(x); print(vv); return(vv); }
  
  # Find the model and iterate until convergence
  model <- Minimzer(model,method="both",lower,upper)
  best <- 1.0e+20
  print(model$fitv,digits=10)
  while (abs(best-model$fitv) > 1.0e-10)
   {
    print("looping")
    print(best)
    best <-model$fitv
    model <- Minimzer(model,method="both",lower,upper)
    print("Objective fn:")
    print(model$fitv,digits=10)
    cat("Difference",best-model$fitv,"\n")
   }
  model <- Minimzer(model,method="both",lower,upper, verbose = verbose)
  print( model$gr(model$env$last.par.best))
  print(model$env$last.par.best)

  # Save tesults
  SaveAll <- NULL
  SaveAll$model <- model
  SaveAll$data <- data
  SaveAll$parameters <- parameters
  SaveAll$report <- model$report()
  save(SaveAll,file=SaveFile)
  rep   <- TMB::sdreport(model)
  if (verbose) print(summary(rep))
  SaveAll$gradient <- rep$gradient.fixed
  SaveAll$sdreport <- rep
  save(SaveAll,file=SaveFile)
  return(model)
 }

# ==============================================================================================================

CreateData <- function(DataFile="data.dat",NDataSet=1, verbose = FALSE)
{
  MatchTable<-function(Table,Char1=NULL,Char2=NULL,Char3=NULL,Char4=NULL,Char5=NULL)
  {
    ii <- rep(T,length(Table[,1]))
    if (!is.null(Char1)) ii <- ii & (Table[,1]==Char1)
    if (!is.null(Char2)) ii <- ii & (Table[,2]==Char2)
    if (!is.null(Char3)) ii <- ii & (Table[,3]==Char3)
    if (!is.null(Char4)) ii <- ii & (Table[,4]==Char4)
    if (!is.null(Char5)) ii <- ii & (Table[,5]==Char5)
    ii <- seq(1:length(Table[,1]))[ii]
    return(ii)
  }

  # Read in the data file
  Data <- read.table(DataFile,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  
  # Extract the minimum and maximum ages
  Index <- MatchTable(Data,Char1="Range_of_ages")
  MinAge <- as.numeric(Data[Index+1,1]); MaxAge <-as.numeric(Data[Index+1,2])   

  # Details of the readers
  IndexVals <- rep(0,NDataSet); Npnt <- rep(0,NDataSet); NReaders <-rep(0,NDataSet)
  MinusA <- rep(0,NDataSet); PlusA <- rep(0,NDataSet); RefAge <- rep(0,NDataSet)
  MaxReader <- 1
  for (Idataset in 1:NDataSet)
   {
    SearchTerm <- paste("Data_set_",Idataset,sep="")
    IndexVals[Idataset] <- MatchTable(Data,Char1=SearchTerm)
    Npnt[Idataset] <- as.numeric(Data[IndexVals[Idataset]+1,1])  
    NReaders[Idataset] <- as.numeric(Data[IndexVals[Idataset]+2,1])  
    MinusA[Idataset] <- as.numeric(Data[IndexVals[Idataset]+3,1])  
    PlusA[Idataset] <- as.numeric(Data[IndexVals[Idataset]+3,2])  
    RefAge[Idataset] <- as.numeric(Data[IndexVals[Idataset]+3,3])  
    Readers <- as.numeric(Data[IndexVals[Idataset]+4,1:NReaders[Idataset]])  
    MaxReader <- max(MaxReader,Readers)
    if (verbose) cat("readers",Readers,"\n")
   }  
  ReadPnt <- matrix(0,nrow=NDataSet,ncol=MaxReader)
  for (Idataset in 1:NDataSet)
    ReadPnt[Idataset,1:NReaders[Idataset]]  <- as.numeric(Data[IndexVals[Idataset]+4,1:NReaders[Idataset]])  

  # Now extract the data
  TheData <- array(-999,dim=c(NDataSet,max(Npnt),MaxReader+1))
  Ipnt <- 0
  for (Idataset in 1:NDataSet)
   { 
    for (Iline in 1:Npnt[Idataset])
     TheData[Idataset,Iline,1:(NReaders[Idataset]+1)] <- as.numeric(Data[(IndexVals[Idataset]+4+Iline),1:(NReaders[Idataset]+1)])
    if (verbose) cat("Last line of data set",Idataset,"is",TheData[Idataset,Npnt[Idataset],1:(NReaders[Idataset]+1)],"\n")
   }

  # Do checks on the data set
  MaAge <- -1; MiAge <- 1000
  NegVals <- 0;
  for (IDataS in 1:NDataSet)
   for (Ipnt in 1:Npnt[IDataS])
    for (Ireader in 1:NReaders[IDataS])
     if (TheData[IDataS,Ipnt,Ireader+1] >= 0)
      {
       if (TheData[IDataS,Ipnt,Ireader+1] > MaAge) MaAge = TheData[IDataS,Ipnt,Ireader+1];
       if (TheData[IDataS,Ipnt,Ireader+1] < MiAge) MiAge = TheData[IDataS,Ipnt,Ireader+1];
      }
    else
     NegVals = 1;
    if (NegVals == 1) cat("WARNING - there are some missing data; the effective sample size calculation may be dubious,\n\n")

  # Create a tabular summary of the data
  write("Structure of the data set",EchoFile,append=T)
  write("Data set # Entries Reader boolean",EchoFile,append=T)
  ReaderStruc <- matrix(0,nrow=1000,ncol=MaxReader+2)
  Presense <- rep(0,MaxReader)
  ReaderSumm <- matrix(0,nrow=NDataSet,ncol=3);
  NrowStruc <- 0;
  for (IDataS in 1:NDataSet)
   for (II in 1:Npnt[IDataS])
    {
     Presense <- rep(0,MaxReader)
     for (Ireader in 1:NReaders[IDataS])
      if (TheData[IDataS,II,Ireader+1] >= 0) Presense[Ireader] = ReadPnt[IDataS,Ireader]
     Ifound <- 0;
     if (NrowStruc>0)
      for (JJ in 1:NrowStruc)
       if (ReaderStruc[JJ,1] == IDataS)
        {
         Jfound = 1;
         for (Ireader in 1:NReaders[IDataS])
          if (Presense[Ireader] != ReaderStruc[JJ,Ireader+2]) Jfound = 0;
         if (Jfound==1) Ifound = JJ;
        }
     if (Ifound==0)
      {
       NrowStruc <- NrowStruc + 1;
       ReaderStruc[NrowStruc,1] = IDataS;
       #ReaderStruc[NrowStruc,2] = 1;
       ReaderStruc[NrowStruc,2] = TheData[IDataS,II,1];
       ReaderSumm[IDataS,3] <- ReaderSumm[IDataS,3] + 1;
       for (Ireader in 1:NReaders[IDataS]) ReaderStruc[NrowStruc,Ireader+2] = Presense[Ireader];
      }
     else
      {
       #ReaderStruc[Ifound,2] <- ReaderStruc[Ifound,2] + 1;
       ReaderStruc[Ifound,2] <- ReaderStruc[Ifound,2] + TheData[IDataS,II,1];
      }
   }
  cat("Number of rows in NrowStruc",NrowStruc,"\n")
  print(NrowStruc)
  ReaderStruc <- matrix(ReaderStruc[1:NrowStruc,],nrow=NrowStruc,ncol=length(ReaderStruc[1,]))
  print("ReaderStruc")
  print(ReaderStruc)
  for (II in 1:NrowStruc) write(ReaderStruc[II,],EchoFile,append=T,ncol=MaxReader+2)
  print("ReaderSumm")
  print(ReaderSumm)

  # Create a reader summary;  ReaderSumm[IDataS,2] is the maximum number of readers for given combination of readers
  for (II in 1:NrowStruc)
   {
    IDataS <- ReaderStruc[II,1];
    ReaderSumm[IDataS,1] = IDataS;
    MaxReaderOpt <- 0;
    for (Ireader in 1:NReaders[IDataS]) if (ReaderStruc[NrowStruc,Ireader+2] >0) MaxReaderOpt <- MaxReaderOpt + 1;
    if (MaxReaderOpt > ReaderSumm[IDataS,2]) ReaderSumm[IDataS,2] = MaxReaderOpt;
   }
  write("ReaderSumm",EchoFile,append=T)
  write(t(ReaderSumm),EchoFile,append=T,ncol=3)
  print("ReaderSumm")
  print(ReaderSumm)

  # Outputs to screen
  write(paste("Number of reads by data set:       ",NReaders),EchoFile,append=T)
  write(paste("Minimum and Maximum Ages:          ",MiAge," ",MaAge),EchoFile,append=T)
  write("",EchoFile,append=T)

  # Record total sample size and specify effective Ns
  EffNOpt <- rep(0,NDataSet)
  TotalN <- rep(0,NDataSet); EffN <- rep(0,NDataSet);
  for (IDataSet in 1:NDataSet)
   {
    TotalN[IDataSet] <- 0;
    for (Ipnt in 1:Npnt[IDataSet]) TotalN[IDataSet] <- TotalN[IDataSet] + TheData[IDataSet,Ipnt,1];
   }
  for (IDataSet in 1:NDataSet)
   if (EffNOpt[IDataSet] <= 0) 
    EffN[IDataSet] <- TotalN[IDataSet] 
  else 
    EffN[IDataSet] <- EffNOpt[IDataSet];

  # Check for duplicates and condense as needed
  OneProblem <- 0;
  for (IDataSet in 1:NDataSet)
   {
    Problem <- 0;
    for (II  in 2:Npnt[IDataSet])
     {
      for (JJ in 1:(II-1))
       if (TheData[IDataSet,JJ,1] > 0)
        {
         Ifound <- 0;
          for (Iread  in 1:NReaders[IDataSet]) if (TheData[IDataSet,JJ,Iread+1] != TheData[IDataSet,II,Iread+1]) Ifound = 1;
          if (Ifound == 0)
           {
            cat("Warning: Lines ",II," and ",JJ," have the same ages\n") 
            TheData[IDataSet,JJ,1] <- TheData[IDataSet,JJ,1] + TheData[IDataSet,II,1];
            TheData[IDataSet,II,1] <- -1;
            Problem <- 1;
           }
        } # JJ
     } 
    if (Problem == 1)
     {
      cat("Duplicate entries found for data set ",IDataSet,"; corrected data set in Echo.File\n")
      cat("Duplicate entries found for data set ", IDataSet, "; corrected data set follows\n")
      NLineOut <- 0;
      for (II in 1:Npnt[IDataSet])
       if (TheData[IDataSet,II,1] > 0) { NLineOut <- NLineOut+1; write(TheData[IDataSet,II,],EchoFile,append=T) } 
      write(paste("New lines ", NLineOut),EchoFile,append=T)
      OneProblem <- 1; 
     }
   }
  if (OneProblem == 1) AA; 
  
  ## Counter for storage
  MaxCells <- 1;
  for (IDataS in 1:NDataSet)
   {
    Ncells <- 1; for (Ireader in 1:ReaderSumm[IDataS,2]) { Ncells <- Ncells * (MaxAge+1); }
    if (Ncells > MaxCells) MaxCells <- Ncells;
   }
  write(paste("total cells ",MaxCells),EchoFile,append=T)

  Outs <- NULL
  Outs$MinAge <- MinAge
  Outs$MaxAge <- MaxAge
  Outs$NDataSet <- NDataSet
  Outs$MaxReader <- MaxReader
  Outs$TheData <- TheData
  Outs$Npnt <- Npnt
  Outs$ReadPnt <- ReadPnt
  Outs$NReaders <- NReaders
  Outs$MinusA <- MinusA
  Outs$PlusA <- PlusA
  Outs$RefAge <- RefAge
  Outs$ReaderSumm <- ReaderSumm
  Outs$ReaderStruc <- ReaderStruc
  Outs$MaxCells <- MaxCells
  Outs$TotalN <- TotalN
  Outs$EffN <- EffN
  if (verbose) print(str(Outs))
  return(Outs)
}
  
# ==============================================================================================================

CreateSpecs <- function(SpecsFile="data.spc",DataSpecs, verbose = FALSE)
{
  MatchTable<-function(Table,Char1=NULL,Char2=NULL,Char3=NULL,Char4=NULL,Char5=NULL)
  {
    ii <- rep(T,length(Table[,1]))
    if (!is.null(Char1)) ii <- ii & (Table[,1]==Char1)
    if (!is.null(Char2)) ii <- ii & (Table[,2]==Char2)
    if (!is.null(Char3)) ii <- ii & (Table[,3]==Char3)
    if (!is.null(Char4)) ii <- ii & (Table[,4]==Char4)
    if (!is.null(Char5)) ii <- ii & (Table[,5]==Char5)
    ii <- seq(1:length(Table[,1]))[ii]
    return(ii)
  }

  Specs <- read.table(SpecsFile,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:10)

  # Read in the initial values  
  ModelSpecs <-vector(mode="list",length=DataSpecs$MaxReader)
  IndexA <- MatchTable(Specs,Char1="#",Char2="reader")+1
  NumBias <- 0; NumSigma <- 0
  for (Ireader in 1:DataSpecs$MaxReader)
   {
    DefaultList <-list(BiasOpt=NULL,SigOpt=NULL,BiasPar=NULL,BiasLow=NULL,BiasHi=NULL,BiasUsed=NULL,SigmaPar=NULL,SigmaLow=NULL,SigmaHi=NULL,SigmaUsed=NULL)  
    DefaultList$BiasOpt <- as.numeric(Specs[IndexA,2])    
    DefaultList$SigOpt <- as.numeric(Specs[IndexA,3])    
    
    # Check for valud bias options
    if (DefaultList$BiasOpt >= 0 & !DefaultList$BiasOpt %in% c(0,1,2)) { cat("Error specifying bias option for reader ",Ireader,"; Bias option ",DefaultList$BiasOpt," is not implemented- stopping\n"); AA }
    if (DefaultList$SigOpt >= 0 & !DefaultList$SigOpt %in% c(1:8)) { cat("Error specifying variance option for reader ",Ireader,"; Variance option ",DefaultList$SigOpt," is not implemented- stopping\n"); AA }
    
    IndexA <- IndexA + 1
    ModelSpecs[[Ireader]] <- DefaultList
   }  

  # Parameters defining Sigmas (Spline)
  Index <- MatchTable(Specs,Char1="#",Char2="Spline",Char3="specifications")
  xvals <- matrix(0,nrow=DataSpecs$MaxReader,ncol=100);
  nknots <- rep(0,DataSpecs$MaxReader)
  for (Ireader in 1:DataSpecs$MaxReader)
   if (ModelSpecs[[Ireader]]$SigOpt == 5)
    {
      nknots[Ireader] <- as.numeric(Specs[Index+1,1])    
      for (IDcnt in 1:nknots[Ireader]) xvals[Ireader,IDcnt] <- as.numeric(Specs[Index+2,IDcnt])
      Index <- Index + 2
      cat("Spline used to define SD for reader ",Ireader,"\n")
      cat("Selected knots are located at "); for (IDcnt in 1:nknots[Ireader]) cat(xvals[Ireader,IDcnt], " "); cat("\n")
    } 
    
  Index <- MatchTable(Specs,Char1="#",Char2="Linear",Char3="specifications")
  # linear model specifications  
  xvalsL <- matrix(0,nrow=DataSpecs$MaxReader,ncol=100);
  nknotsL <- rep(0,DataSpecs$MaxReader)
  for (Ireader in 1:DataSpecs$MaxReader)
    if (ModelSpecs[[Ireader]]$SigOpt == 6)
    {
      nknotsL[Ireader] <- as.numeric(Specs[Index+1,1])    
      for (IDcnt in 1:nknotsL[Ireader]) xvalsL[Ireader,IDcnt] <- as.numeric(Specs[Index+2,IDcnt])
      Index <- Index + 2
      if (xvalsL[Ireader,1] != DataSpecs$MinAge) { print("First age must be 1"); AA }
      if (xvalsL[Ireader,nknotsL[Ireader]] != DataSpecs$MaxAge) { print("last age must be MaxAge"); AA; }
      cat("Linear interpolation used to define SD for reader ",Ireader,"\n")
      cat("Selected knots are located at "); for (IDcnt in 1:nknotsL[Ireader]) cat(xvalsL[Ireader,IDcnt], " "); cat("\n")
    } 

  # Read in the initial values  
  IndexB <- MatchTable(Specs,Char1="Bias_Pars")+1
  IndexC <- MatchTable(Specs,Char1="Sigma_Pars")+1
  NumBias <- 0; NumSigma <- 0
  for (Ireader in 1:DataSpecs$MaxReader)
   {
    DefaultList <-list(BiasOpt=NULL,SigOpt=NULL,BiasPar=NULL,BiasLow=NULL,BiasHi=NULL,BiasUsed=NULL,SigmaPar=NULL,SigmaLow=NULL,SigmaHi=NULL,SigmaUsed=NULL)  
    DefaultList$BiasOpt <- ModelSpecs[[Ireader]]$BiasOpt
    DefaultList$SigOpt <- ModelSpecs[[Ireader]]$SigOpt
    Nbias <- 0
    if (DefaultList$BiasOpt==1) Nbias <- 1
    if (DefaultList$BiasOpt==2) Nbias <- 3
    if (Nbias>0)
     for (Ibias in 1:Nbias)
      {
       DefaultList$BiasPar <- c(DefaultList$BiasPar,as.numeric(Specs[IndexB,3]))
       DefaultList$BiasLow <- c(DefaultList$BiasLow,as.numeric(Specs[IndexB,1]))
       DefaultList$BiasHi <- c(DefaultList$BiasHi,as.numeric(Specs[IndexB,2]))
       DefaultList$BiasUsed <- c(DefaultList$BiasUsed,as.numeric(Specs[IndexB,4]))
       IndexB <- IndexB + 1
      }  
    NumBias <- NumBias + Nbias
    Nsigma <- 0
    if (DefaultList$SigOpt==1) Nsigma <- 1
    if (DefaultList$SigOpt==2) Nsigma <- 3
    if (DefaultList$SigOpt==3) Nsigma <- 3
    if (DefaultList$SigOpt==5) Nsigma <- nknots[Ireader];
    if (DefaultList$SigOpt==6) Nsigma <- nknotsL[Ireader];
    if (DefaultList$SigOpt==7) Nsigma <- 2
    if (DefaultList$SigOpt==8) Nsigma <- 2
    if (Nsigma>0)
     for (Isigma in 1:Nsigma)
      {
       DefaultList$SigmaPar <- c(DefaultList$SigmaPar,as.numeric(Specs[IndexC,3]))
       DefaultList$SigmaLow <- c(DefaultList$SigmaLow,as.numeric(Specs[IndexC,1]))
       DefaultList$SigmaHi <- c(DefaultList$SigmaHi,as.numeric(Specs[IndexC,2]))
       DefaultList$SigmaUsed <- c(DefaultList$SigmaUsed,as.numeric(Specs[IndexC,4]))
       IndexC <- IndexC + 1
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
  if (verbose) print(str(Outs))
  return(Outs)

}
# ==============================================================================================================

#' Plot output
#'
#' Plots age comparisons and results from the fitted Ageing Error model
#'
#' @param Data Input data matrix
#' @param MaxAge Maximum estimated age
#' @param SaveFile Directory for fitted model
#' @param PlotType Type of saved plots, i.e. PDF or PNG
#' @param ReaderNames Vector with names of each reader, defaults to
#' "Reader 1", "Reader 2", etc.
#' @param subplot Vector of which plots to create.
#' @param dots Additional arguments passed to the
#' \code{\link{ageing_comparison}} function.
#' @return Returns AIC, AICc, and BIC for fitted model.
#'
#' @references Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australias southern and eastern scalefish
#' and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
#'
#' @author James T. Thorson, Ian G. Taylor
#'
#' @export
#'
PlotOutputFn <-
  function(Data, IDataSet, MaxAge, Report, PlotType = "PNG", SaveFile=getwd(), subplot=1:3, Nparameters=0,LogLike=0,
           ReaderNames = NULL, Species=NULL, SaveDir="", verbose = FALSE, ...)
  {
    SaveFile <- paste(SaveFile,"/",SaveDir,sep="")

    # Interpret inputs
    Nreaders <- ncol(Data)-1
    Ages <- Nages <- MaxAge+1

    # Reader names
    if( is.null(ReaderNames) ){
      ReaderNames <- paste("Reader", 1:Nreaders)
    }
    
    # Age-reading error matrices: dimensions are Reader, TrueAge, EstAge
    MisclassArray <- array(NA, dim = c(Nreaders, Ages, Ages),
                           dimnames = list(paste("Reader", 1:Nreaders),
                                           paste("TrueAge", 0:MaxAge),
                                           paste("EstAge", 0:MaxAge)))
    for(i in 1:Nreaders){
      MisclassArray[i,,] <-  Report$AgeErrOut[i,,]
    }

    # Estimated age-structure
    AgeStruct <- cbind(0:MaxAge,t(Report$Aprob))

    # Reader CV, SD and Bias
    Temp <- matrix(0,nrow=5*Nages*Nreaders,ncol=5)
    for (Ireader in 1:Nreaders)
     {
      yrange <- (Ireader-1)*Nages
      Temp[yrange+1:Nages,1]  <- Ireader
      Temp[yrange+1:Nages,2]  <- 0:MaxAge
      Temp[yrange+1:Nages,3]  <- Report$TheSD[Ireader,]/c(1,1:MaxAge)
      Temp[yrange+1:Nages,4]  <- Report$TheSD[Ireader,]
      Temp[yrange+1:Nages,5]  <- Report$TheBias[Ireader,]
     }  
    Temp <- t(Temp)

    ErrorAndBiasArray <- array(as.numeric(Temp),
                               dim = c(5, Nages, Nreaders),
                               dimnames = list(c("Reader", "True_Age",
                                                 "CV", "SD", "Expected_age"),
                                               paste("Age", 0:MaxAge),
                                               paste("Reader", 1:Nreaders)))
    if (verbose) print(str(ErrorAndBiasArray))
    # Estimate unobserved age for each otolith
    # This is done by assigning each otolith to the age which has
    # maximum posterior probability (i.e. the conditional mode,
    # as is typically done for random effects)
    AgeProbs <- array(NA,
                      dim = c(nrow(Data), Ages),
                      dimnames = list(paste("Otolith", 1:nrow(Data)),
                                      paste("TrueAge", 0:MaxAge)))
    OtI <- AgeI <- ReadI <- 1
    for(OtI in 1:nrow(Data)){
      for(AgeI in 1:Ages){
        AgeProbs[OtI, AgeI] <- 1
        for(ReadI in 1:Nreaders){
          if(Data[OtI, ReadI+1] != -999){
            AgeRead <- Data[OtI, ReadI+1]
            AgeProbs[OtI, AgeI] <- AgeProbs[OtI, AgeI] *
              (MisclassArray[ReadI, AgeI, AgeRead+1])^Data[OtI, 1]
          } # end check for value other than -999
        } # end loop over readers
      } # end loop over ages
    } # end loop over rows of data
    
    # Remove MaxAge before calculating "TrueAge"
    # because the MaxAge is a plus-group, and ends up with
    # maximum probability for most ages in the upper tail
    # ANDRE - Found another issue - should be age-1 because the first age is zero
    TrueAge <- apply(AgeProbs, MARGIN = 1,
                     FUN = function(Vec){order(Vec[-length(Vec)],
                                               decreasing = TRUE)[1]})-1  
    
    DataExpanded <- Data[rep(1:nrow(Data), Data[, 1]), -1]
    DataExpanded[DataExpanded == -999] <- NA
    
    ####################################################################
    # Plot comparison of data for each pair of readers
    

    if(1 %in% subplot){
      if(PlotType == "PDF"){
        pdf(paste0(SaveDir,"/",Species,"-Data set",IDataSet,ReaderNames[ireader],
                   " vs ", ReaderNames[jreader], ".pdf",sep=""),
            width = 6, height = 6)
      }
      
      # make plots of input data for each reader pair
      for(ireader in 1:(Nreaders-1)){
        for(jreader in (ireader+1):Nreaders){
          ageing_comparison(xvec = DataExpanded[,ireader],
                            yvec = DataExpanded[,jreader],
                            xlab = ReaderNames[ireader],
                            ylab = ReaderNames[jreader],
                            maxage = max(DataExpanded, na.rm=TRUE),
                            hist=F,
                            png = (PlotType == "PNG"),
                            SaveFile = SaveFile,
                            filename = paste0(Species,"-Data set-",IDataSet," ",ReaderNames[ireader],
                                              " vs ", ReaderNames[jreader], ".png",sep=""),
                            verbose = verbose,
                            ...)
        }
      }
    } # end check for whether subplot 1 was requested

    ####################################################################
    # Plot estimated age structure
    
    if(2 %in% subplot){
      if(PlotType == "PDF"){
        pdf(file.path(SaveFile, paste(Species,"-Data set-",IDataSet,"Estimated vs Observed Age Structure.pdf",sep="")),
            width = 6, height = 6)
      }
      if(PlotType == "PNG"){
        png(file.path(SaveFile, paste(Species,"-Data set-",IDataSet,"Estimated vs Observed Age Structure.png",sep="")),
            width = 6, height = 6, units = "in", res = 200)
      }
      par(mar = c(3, 3, 2, 0), mgp = c(1.5, 0.25, 0),
          tck = -0.02, oma = c(0, 0, 0, 0)+0.1)
      plot(x = AgeStruct[, 1], y = AgeStruct[, 2], type = "s", lwd = 2,
           xlab = "Age", ylab = "Prop", main = "Estimated=Black, Observed=Red")
      hist(as.matrix(DataExpanded),
           add = TRUE, freq = FALSE, breaks = seq(0, MaxAge, by = 1),
           col = rgb(red=1, green=0, blue=0, alpha=0.30))
      dev.off()
    } # end check for whether subplot was requested
    
    ####################################################################
    # Plot true age against different age reads
    
    if(3 %in% subplot){
      Ncol <- ceiling(sqrt(Nreaders))
      Nrow <- ceiling(Nreaders/Ncol)
      if(PlotType == "PDF"){
        pdf(file.path(SaveFile, paste(Species,"-Data set-",IDataSet,"True vs Reads (by reader).pdf",sep="")),
            width = Ncol*3, height = Nrow*3)
      }
      if(PlotType == "PNG"){
        png(file.path(SaveFile, paste(Species,"-Data set-",IDataSet,"True vs Reads (by reader).png",sep="")),
            width = Ncol*3, height = Nrow*3, units = "in", res = 200)
      }
      par(mfrow = c(Nrow, Ncol), mar = c(3, 3, 2, 0), mgp = c(1.5, 0.25, 0),
          tck = -0.02, oma = c(0, 0, 5, 0)+0.1)
      for(ReadI in 1:Nreaders){
        Main <- ReaderNames[ReadI]
        
        # Add 0.5 to match convention in Punt model that otoliths are read
        # half way through year
        Temp <- cbind(TrueAge, Data[, ReadI+1]+0.5)
        # Exclude rows with no read for this reader
        Temp <- Temp[which(Data[, ReadI+1] != -999), ]
        plot(x = Temp[, 1], y = Temp[, 2],
             ylim = c(0, MaxAge), xlim = c(0, MaxAge),
             col = rgb(red=0, green=0, blue=0, alpha=0.2),
             xlab = "Mode predicted age | parameters",
             ylab = "Read age", lwd = 2, main = Main, pch = 21, cex = 0.2)
        lines(x = c(0, MaxAge), y = c(0, MaxAge), lwd = 1, lty = "dashed")
        lines(x = ErrorAndBiasArray['True_Age', , ReadI],
              y = ErrorAndBiasArray['Expected_age', , ReadI],
              type = "l", col = "red", lwd = 1)
        lines(x = ErrorAndBiasArray['True_Age', , ReadI],
              y = ErrorAndBiasArray['SD', , ReadI],
              type = "l", col = "blue", lwd = 1)
        lines(x = ErrorAndBiasArray['True_Age', , ReadI],
              y = ErrorAndBiasArray['Expected_age', , ReadI] +
                2*ErrorAndBiasArray['SD', , ReadI],
              type = "l", col = "red", lwd = 1, lty = "dashed")
        lines(x = ErrorAndBiasArray['True_Age', , ReadI],
              y = ErrorAndBiasArray['Expected_age', , ReadI] -
                2*ErrorAndBiasArray['SD', , ReadI],
              type = "l", col = "red", lwd = 1, lty = "dashed")
      }
      mtext(side = 3, outer = TRUE,
            text = paste0("Reads(dot), Sd(blue), expected_read(red solid line),\n",
                          " and 95% CI for expected_read(red dotted line)"),
            line = 1)
      dev.off()
    } # end check for whether subplot was requested
    
    ## AIC
    Nll <- LogLike
    Df <- Nparameters
    n <- sum(ifelse(Data[, -1] == -999, 0, 1))
    Aic <- 2*Nll + 2*Df
    Aicc <- Aic + 2*Df*(Df+1)/(n-Df-1)
    Bic <- 2*Nll + Df*log(n)
    
    # Write definitions to file
    for(ReadI in 1:Nreaders){
      Main <- ReaderNames[ReadI]
      write.csv( ErrorAndBiasArray[, , ReadI],
                 file = file.path(SaveFile, paste0(Species," SS_format_", Main, ".csv",sep="")))
    }
    
    
    # Return stuff
    ModelSelection <- list("AIC" = Aic,  "AICc" = Aicc,  "BIC" = Bic)
    Output <- list("ModelSelection" = ModelSelection,
                   "ErrorAndBiasArray" = ErrorAndBiasArray)
    return(Output)
  }

# ==============================================================================================================

ProcessResults <- function(Species,SaveDir,CalcEff=F, verbose = FALSE)
{
 
  
 SaveFile <- paste(SaveDir,"/",Species,".lda",sep="")
 print(SaveFile)
 ReportFile <- paste(SaveDir,"/",Species,".rpt",sep="")
 print(ReportFile)
 
 load(SaveFile)
 if (verbose) print(str(SaveAll))
 NDataSet <- SaveAll$data$NDataSet
 MaxReader <- SaveAll$data$MaxReader
 
 write(paste(SaveAll$report$f," ",SaveAll$report$Obj_fun),ReportFile)
 write(paste("Total number of readers:",MaxReader),ReportFile,append=T)
 write(paste("Number of data sets:",NDataSet),ReportFile,append=T)
 write(paste("Bias options by reader:",SaveAll$data$BiasOpt),ReportFile,append=T)
 write(paste("Sigma options by reader:",SaveAll$data$SigOpt),ReportFile,append=T)
 
 write(paste("Total objective function:",SaveAll$report$f),ReportFile,append=T)
 Index <- which(abs(SaveAll$gradient)== max(abs(SaveAll$gradient)))
 write(paste("maximum gradient:",SaveAll$gradient[Index]),ReportFile,append=T)
 
 write(paste("Number of readers: ",MaxReader),ReportFile,append=T)
 write(paste("Range of ages: ", SaveAll$data$MinAge , " - ", SaveAll$data$MaxAge),ReportFile,append=T)
 for (Ireader in 1:NDataSet)
   write(paste("Reader #", Ireader, " Minus/Plus ages: ", SaveAll$data$MinusA[Ireader], " / ", SaveAll$data$PlusA[Ireader]),ReportFile,append=T)
 write(paste("Number of data sets:", NDataSet),ReportFile,append=T)
 write(paste("Number of lines of data per data set: ", SaveAll$data$Npnt),ReportFile,append=T)
 write(paste("Number of data points per data set:", SaveAll$data$TotalN),ReportFile,append=T)
 write("\n",ReportFile,append=T)
 write("# Bias options",ReportFile,append=T)
 write("#   -X: Mirrored with pattrn X",ReportFile,append=T)
 write("#   0: Unbiased [0 parameters]",ReportFile,append=T)
 write("#   1: Linearly proportional to age [1 parameter]",ReportFile,append=T)
 write("#   2: Michaelis-Menten [3 parameters]",ReportFile,append=T)
 write("\n# SD options",ReportFile,append=T)
 write("#   1: Constant CV [1 parameter]",ReportFile,append=T)
 write("#   2: SD Michaelis-Menten [3 parameters]",ReportFile,append=T)
 write("#   3: CV Michaelis-Menten [3 parameters]",ReportFile,append=T)
 write("#   4: Known perfectly [0 parameters]",ReportFile,append=T)
 write("#   5: SD spline function of age [variable parameters]",ReportFile,append=T)
 write("#   6: SD linear piecewise function of age [variable parameters]",ReportFile,append=T)
 write("#   7: SD linear function of age [2 parameters]",ReportFile,append=T)
 write("#   8: CV linear function of age [2 parameters]",ReportFile,append=T)
 write(paste("Reader BiasType SigmaType"),ReportFile,append=T);
 for (Ireader in 1:MaxReader)
   write(paste(Ireader, " ", SaveAll$data$BiasOpt[Ireader], " ", SaveAll$data$SigOpt[Ireader]),ReportFile,append=T)
 write("",ReportFile,append=T)
 
 # Bias and variance 
 write("Reader Age CV SD Expected age",ReportFile,append=T)
 for (Ireader in 1:MaxReader)
  for (Age in 0:SaveAll$data$MaxAge)
   {
    if (Age > 1)
      CVV <- SaveAll$report$TheSD[Ireader,Age+1]/(Age*1.0)
    else
      CVV <- SaveAll$report$TheSD[Ireader,Age+1]
    write(paste(Ireader," ", Age," ", round(CVV,5)," ",round(SaveAll$report$TheSD[Ireader,Age+1],5)," ", round(SaveAll$report$TheBias[Ireader,Age+1],5)),ReportFile,append=T)
   }
 write("",ReportFile,append=T)
 
 # Estimated age-structure
 write("Estimated age-structure by data set",ReportFile,append=T)
 write("===================================",ReportFile,append=T)
 HeadString <- "Age "; for (IDataSet in 1:NDataSet) HeadString <- paste(HeadString,"Data set#",IDataSet)
 write(HeadString,ReportFile,append=T)
 for (Age in 0:SaveAll$data$MaxAge)
  {
   HeadString <- "Age "; for (IDataSet in 1:NDataSet) HeadString <- paste(HeadString,SaveAll$report$Aprob[IDataSet,Age+1])
   write(HeadString,ReportFile,append=T)
  }  
 write("",ReportFile,append=T)
 
 # Age-reading error matrices
 write("Final age-reading error matrices",ReportFile,append=T)
 for (Ireader in 1:MaxReader) 
  {
   write(paste("Matrix for reader# ",Ireader),ReportFile,append=T)
   write(t(SaveAll$report$AgeErrOut[Ireader,,]),ncol=SaveAll$data$MaxAge+1,ReportFile,append=T)
 }
 
 if (!is.null(SaveAll$sdreport))
  {
   write("\nVariable estimate SD",ReportFile,append=T)
   names <- row.names(summary(SaveAll$sdreport))
   write(t(cbind(names,summary(SaveAll$sdreport))),ReportFile,append=T,ncol=3)
  }

 # find the effective sample sizes
 if (CalcEff==T)
  {
   ReaderSumm <- SaveAll$data$ReaderSumm
   ReaderStruc <- SaveAll$data$ReaderStruc
   write("Compute the effective sample sizes",ReportFile,append=T)
   write("==================================",ReportFile,append=T)
   GroupPointer <- 0
   Ages <- rep(0,MaxReader)
   ProbStore <- rep(0,SaveAll$data$MaxCells)
   for (IDataSet in 1:NDataSet)
    {
     write(paste("Data set: ", IDataSet),ReportFile,append=T)
     write(paste("Data_set Group Group Line Readers Obs Obs_Numbers Pred_Numbers"),ReportFile,append=T)
     Top <- 0; Bot <- 0;
     for (Kgroup in 1:ReaderSumm[IDataSet,3])
      {
       GroupPointer <- GroupPointer + 1
 
       # Find the total number of combinations of ages (MaxAge+1)**number_of_readers
       Ncells <- 1; for (Iread in 1:ReaderSumm[IDataSet,2]) Ncells <- Ncells * (SaveAll$data$MaxAge+1);

       # Move through each possible combination of ages between 0 and MaxAge 
       TotalProb <- 0;
       for (II in 1:Ncells)
        {
       
         # Find the ages for this "cell" 
         Ndiv <- II
         for (Iread in 1:SaveAll$data$Nread[IDataSet]) Ages[Iread] <- -1;
         for (Iread in 1:SaveAll$data$Nread[IDataSet])
          if (ReaderStruc[GroupPointer,Iread+2] > 0)
           {
            DivJ <- 1; 
            if (Iread < SaveAll$data$Nread[IDataSet])
             for (Jread in SaveAll$data$Nread[IDataSet]:(Iread+1)) 
              if (Jread>0 & Jread<=SaveAll$data$Nread[IDataSet])
               if (ReaderStruc[GroupPointer,Jread+2] > 0) DivJ <- DivJ * (SaveAll$data$MaxAge+1);
            DivI <- floor((Ndiv-1)/DivJ)
            Ages[Iread] <- DivI
            Ndiv <- Ndiv - DivI*DivJ
          }

         # Find the probability for this cell, i.e. the probability of an ageing reading of Ages(1)&Ages(2)&...
         if (SaveAll$data$Iperfect[IDataSet] == 0)
          {
           Prob1 = 0;
           for (Age1 in SaveAll$data$MinAge:SaveAll$data$MaxAge)
            {
             # Prior probability * product over readers
             Prob2 <- SaveAll$report$Aprob[IDataSet,Age1+1]
             for (Ireader in 1:SaveAll$data$Nread[IDataSet])
              { 
               Jreader <- SaveAll$data$ReadPnt[IDataSet,Ireader];
               AgeA = Ages[Ireader]; 
               if (AgeA >= 0) Prob2 = Prob2 * SaveAll$report$AgeErrOut[Jreader,Age1+1,AgeA+1] 
              } 
             Prob1 <- Prob1 + Prob2; 
            }
           ProbStore[II] <- Prob1; 
           TotalProb <- TotalProb + Prob1;   
          }
         else
          {
           # Product over readers
           Age1 = Ages[1];
           Prob1 <- 1;
           for (Ireader in 2:SaveAll$data$Nread[IDataSet])
            { 
             Jreader <- SaveAll$data$ReadPnt[IDataSet,Ireader];
             AgeA <- Ages[Ireader]
             if (Age1 >=0 & AgeA >= 0) Prob1 = Prob1 * SaveAll$report$AgeErrOut[Ireader,Age1+1,AgeA+1]; 
            }
           ProbStore[II] <- Prob1; 
           TotalProb <- TotalProb + Prob1;   
          }
        } # for (II in 1:Ncells)

       # Now compute
       writeout <- NULL
       for (II in 1:Ncells)
        {
       
         # Find the ages for this "cell" 
         Ndiv <- II;
         for (Iread in 1:SaveAll$data$Nread[IDataSet]) Ages[Iread] <- -1;
         for (Iread in 1:SaveAll$data$Nread[IDataSet])
           if (ReaderStruc[GroupPointer,Iread+2] > 0)
            {
             DivJ <- 1; 
             if (Iread < SaveAll$data$Nread[IDataSet])
             for (Jread in SaveAll$data$Nread[IDataSet]:(Iread+1)) 
              if (Jread>0 & Jread<=SaveAll$data$Nread[IDataSet])
               if (ReaderStruc[GroupPointer,Jread+2] > 0) DivJ <- DivJ * (SaveAll$data$MaxAge+1);
             DivI <- floor((Ndiv-1)/DivJ)
             Ages[Iread] <- DivI
             Ndiv <- Ndiv - DivI*DivJ
           }

         # Check for a match
         Jfound <- 0; Pobs = 0;
         for (JJ in 1:SaveAll$data$Npnt[IDataSet])
          {
           # Set Ifound to 1 if the ages don't match properly
           Ifound <- 0;
           for (Iread in 1:SaveAll$data$Nread[IDataSet])
             if (Ages[Iread] != SaveAll$data$TheData[IDataSet,JJ,Iread+1] & SaveAll$data$TheData[IDataSet,JJ,Iread+1] >=0) Ifound <- 1;
           # We have a match so store the results 
           if (Ifound == 0)
            {  Pobs = SaveAll$data$TheData[IDataSet,JJ,1]/ReaderStruc[GroupPointer,2]; Jfound <- 1; KK <- JJ; }
          }

         # Find the probability for this cell 
         Pest = ProbStore[II] / TotalProb+1.0e-190; 
         if (Jfound == 1) 
          {
           HeadString <- paste("Data Point: ", IDataSet, " ", GroupPointer, " ", Kgroup, " ", KK, " ")
           for (Iread in 1:SaveAll$data$Nread[IDataSet]) HeadString <- paste(HeadString,Ages[Iread], " ")
           HeadString <- paste(HeadString,SaveAll$data$TheData[IDataSet,KK,1], " ", round(Pobs*SaveAll$data$TotalN[IDataSet],7), " ", round(Pest*SaveAll$data$TotalN[IDataSet],12), " ")
           writeout <- rbind(writeout,HeadString)   
          } 

         # Compute the effective sample size muliplier
         Top <- Top + Pest*(1-Pest)/ReaderStruc[GroupPointer,2];
         Bot <- Bot + (Pest-Pobs)^2;
       }
      write(writeout,ReportFile,append=T)
     } # Kgroup
    EffPred = Top/Bot*SaveAll$data$TotalN[IDataSet];
    write("Data_Set Predicted_EFF Assumed_Eff Sample_size",ReportFile,append=T)
    write(paste(IDataSet, " ", EffPred, " ", SaveAll$data$EffN[IDataSet], " ", SaveAll$data$TotalN[IDataSet]),ReportFile,append=T)
    write("",ReportFile,append=T)
   }

  } # IDataSet
 Npars <- length(SaveAll$sdreport$par.fixed)
 
 # find the effective sample sizes
 for (IDataSet in 1:NDataSet)
  {
   Data <- SaveAll$data$TheData[IDataSet,1:SaveAll$data$Npnt[IDataSet],]
   Output <- PlotOutputFn(Data, IDataSet,SaveAll$data$MaxAge, SaveAll$report, Nparameters=Npars,LogLike=SaveAll$report$Obj_fun,PlotType = "PNG", subplot=1:3, ReaderNames = NULL, Species=Species,SaveDir=SaveDir)
  }   
 #print(str(Output))
 return(Output)
 
}
