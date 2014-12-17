
RunFn <-
function(Data, SigOpt, KnotAges, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize=0, Intern=TRUE, AdmbFile=NULL, JustWrite=FALSE, CallType="system", ExtraArgs=" -est"){

  # Copy ADMB file 
  if(!is.null(AdmbFile)) file.copy(from=paste(AdmbFile,"agemat.exe",sep=""), to=paste(SaveFile,"agemat.exe",sep=""), overwrite=TRUE)
  
  # Check for errors
  Nreaders = ncol(Data)-1
  for(ReaderI in 1:Nreaders){
    if( (SigOpt[ReaderI]==5 | SigOpt[ReaderI]==6) & is.na(KnotAges[[ReaderI]][1]) ) stop("Must specify KnotAges for any reader with SigOpt 5 or 6")
  } 
  
  # Check for specification errors
  for(ReaderI in 1:Nreaders){
    if( (SigOpt[ReaderI]<0 & SigOpt[ReaderI]<=(-ReaderI)) | (BiasOpt[ReaderI]<0 & BiasOpt[ReaderI]<=(-ReaderI)) ) stop("Mirrored readers must mirror a lower numbered reader")
  }
  
  # Write DAT file
  write(c("# Maximum number of readers",Nreaders),file=paste(SaveFile,"agemat.dat",sep=""))
    write(c("# Number of data sets",NDataSets),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Number of points per data set",nrow(Data)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Readers per data set",ncol(Data)-1),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write("# Which readers per data set",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(1:(ncol(Data)-1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write(c("# Minimum age",MinAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Maximum age",MaxAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Reference age",RefAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Minus groups",MinusAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Plus groups",PlusAge),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write("# Option for bias",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(BiasOpt),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write("# Option for standard deviation",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(SigOpt),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE,col.names=FALSE)
    write(c("# Option for effective sample size",EffSampleSize),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write(c("# Use Par File (1=Yes)",0),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE) 
  # Write knots related to splines
    write("\n# Number and location of knots for any splines",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(ReaderI in 1:Nreaders){
      if(SigOpt[ReaderI]==5 | SigOpt[ReaderI]==6) write.table( rMx(c(length(KnotAges[[ReaderI]]),KnotAges[[ReaderI]])), file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE) 
    }
  # Write initial values  
    # Bias 
    write("\n# Min, Max, Init, Phase for Bias",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(BiasI in 1:Nreaders){
      # No bias
      if(BiasOpt[BiasI]<=0){}
      # Linear bias
      if(BiasOpt[BiasI]==1) write.table(rMx(c(0.001,3,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      # Curvilinear bias = 0.5+Par1 + (Par3-Par1)/(1.0-mfexp(-Par2*(float(MaxAge)-1)))*(1.0-mfexp(-Par2*(float(Age1)-1)))
      # Starting value must be non-zero        
      if(BiasOpt[BiasI]==2){
        write.table(rMx(c(0.001,10,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.01,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,MaxAge*2,MaxAge,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
    }
    # Sigma
    write("\n# Min, Max, Init, Phase for Sigma",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(SigI in 1:Nreaders){
      # No error
      if(SigOpt[SigI]<=0){}
      # Linear CV
      if(SigOpt[SigI]==1) write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      # Curvilinear SD = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
      # Starting value must be non-zero
      if(SigOpt[SigI]==2){
        write.table(rMx(c(0.001,100,1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.01,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,100,10,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
      # Curvilinear CV
      # Curvilinear CV = Par1 + (Par3-Par1)/(1.0-exp(-Par2*(100-1)))*(1.0-exp(-Par2*(1:100)))
      # Starting value must be non-zero
      if(SigOpt[SigI]==3){
        write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(-10,1,0.01,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        write.table(rMx(c(0.001,3,0.1,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      }
      # Spline with estimated derivative at beginning and end (Params 1-N: knot parameters; N+1 and N+2: derivative at beginning and end)
      if(SigOpt[SigI]==5){
        for(ParI in 1:(2+length(KnotAges[[SigI]]))){
          write.table(rMx(c(-10.0,10.0,1.0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        }
      }
      # Spline with derivative at beginning and end fixed at zero (Params 1-N: knot parameters)
      if(SigOpt[SigI]==6){
        for(ParI in 1:length(KnotAges[[SigI]])){
          write.table(rMx(c(-10.0,10.0,1.0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
        }
      }
    }
    # Probs (i.e. age-composition probability relative to reference age)
    write("\n# Min, Max, Phase for Probs",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(rMx(c(-20,20,2)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)
    # Slopes
    write("\n# Min, Max, Init, Phase for slopes",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    for(DataSetI in 1:NDataSets){
      if(MinAge>1) write.table(rMx(c(-10,0,0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
      write.table(rMx(c(-10,0,0,1)),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,col.names=FALSE,row.names=FALSE)      
    }    
  # Write dataset    
    write("\n# Data set",file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)
    write.table(Data,file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE,row.names=FALSE, col.names=FALSE)
    write(c("# Test number",123456),file=paste(SaveFile,"agemat.dat",sep=""),append=TRUE)

  # Run ADMB file
  if(JustWrite==FALSE){
    setwd(SaveFile)
    if(CallType=="shell") Output = shell( paste0("agemat.exe",ExtraArgs),intern=Intern)   # This may need to have the location pasted onto it depending upon file structure
    if(CallType=="system") Output = system( paste0("agemat.exe",ExtraArgs),intern=Intern)
    #Admb = scan(paste(SaveFile,"agemat.par",sep=""),comment.char="#",quiet=TRUE)
  }
  #return(Output)
}

