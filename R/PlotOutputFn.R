
PlotOutputFn <-
function(Data, MaxAge, SaveFile, PlotType="PDF", ReaderNames=NULL){

  # Interpret inputs
  Nreaders = ncol(Data)-1
  Ages = Nages = MaxAge+1
  
  # Read REP file
  Rep = scan(paste(SaveFile,"agemat.rep",sep=""),comment.char="%", what="character", quiet=TRUE)
    
  # Read Misclassification rates
  Grep = grep("reader#", Rep)
    MisclassArray = array(NA, dim=c(Nreaders,Ages,Ages), dimnames=list(paste("Reader",1:Nreaders),paste("TrueAge",0:MaxAge),paste("EstAge",0:MaxAge)))   # Reader, TrueAge, EstAge
    for(i in 1:Nreaders){
      MisclassArray[i,,] = matrix(as.numeric(Rep[Grep[i]+1+1:(Ages^2)]),ncol=Ages,byrow=TRUE)
    }
  
  # Input estimated age-structure
  Grep = grep("age-structure", Rep)
    AgeStruct = matrix(as.numeric(Rep[Grep[1]+7+1:(2*Ages)]),ncol=2,byrow=TRUE)
    
  # Input reader error and bias  
  Grep = grep("age",Rep)[3]
    Temp = Rep[Grep+1:(5*Nages*Nreaders)]
    ErrorAndBiasArray = array(as.numeric(Temp), dim=c(5,Nages,Nreaders), dimnames=list(c("Reader","True_Age","CV","SD","Expected_age"),paste("Age",0:MaxAge),paste("Reader",1:Nreaders)))
  
  # Estimate unobserved age for each otolith
    # This is done by assigning each otolith to the age which has maximum posterior probability (i.e. the conditional mode, as is typically done for random effects)
  AgeProbs = array(NA, dim=c(nrow(Data),Ages), dimnames=list(paste("Otolith",1:nrow(Data)),paste("TrueAge",0:MaxAge)))
    OtI = AgeI = ReadI = 1
  for(OtI in 1:nrow(Data)){
  for(AgeI in 1:Ages){
    AgeProbs[OtI,AgeI] = AgeStruct[AgeI,2]
    for(ReadI in 1:Nreaders){
      if(Data[OtI,ReadI+1]!=-999){
        AgeRead = Data[OtI,ReadI+1]
        AgeProbs[OtI,AgeI] = AgeStruct[AgeI,2] * (MisclassArray[ReadI,AgeI,AgeRead+1])^Data[OtI,1]
      }
    }
  }}
  # Remove MaxAge before calculating "TrueAge" because the MaxAge is a plus-group, and ends up with maximum probability for most ages in the upper tail
  TrueAge = apply(AgeProbs, MARGIN=1, FUN=function(Vec){order(Vec[-length(Vec)],decreasing=TRUE)[1]})
  
  # Plot estimated age structure
  if(PlotType=="PDF") pdf(paste(SaveFile,"Estimated vs Observed Age Structure.pdf",sep=""),width=6,height=6)
  if(PlotType=="PNG") png(paste(SaveFile,"Estimated vs Observed Age Structure.png",sep=""),width=6,height=6,units="in",res=200)
    par(mar=c(3,3,2,0),mgp=c(1.5,0.25,0),tck=-0.02,oma=c(0,0,0,0)+0.1)
    plot(x=AgeStruct[,1],y=AgeStruct[,2],type="s",lwd=2,xlab="Age",ylab="Prop",main="Estimated=Black, Observed=Red")
    DataExpanded = Data[rep(1:nrow(Data),Data[,1]),-1]
    hist(ifelse(DataExpanded==-999,NA,DataExpanded),add=TRUE,freq=FALSE,breaks=seq(0,MaxAge,by=1),col=rgb(red=1,green=0,blue=0,alpha=0.30))
  dev.off()
    
  # Plot true age against different age reads
  Ncol=ceiling(sqrt(Nreaders))
    Nrow=ceiling(Nreaders/Ncol)
  if(PlotType=="PDF") pdf(paste(SaveFile,"True vs Reads (by reader).pdf",sep=""),width=Ncol*3,height=Nrow*3)
  if(PlotType=="PNG") png(paste(SaveFile,"True vs Reads (by reader).png",sep=""),width=Ncol*3,height=Nrow*3,units="in",res=200)
    par(mfrow=c(Nrow,Ncol),mar=c(3,3,2,0),mgp=c(1.5,0.25,0),tck=-0.02,oma=c(0,0,5,0)+0.1)
    for(ReadI in 1:Nreaders){
      if( is.null(ReaderNames) ){
        Main = paste("Reader",ReadI)
      }else{
        Main = ReaderNames[ReadI]
      }
      Temp = cbind(TrueAge, Data[,ReadI+1]+0.5)   # Add 0.5 to match convention in Punt model that otoliths are read half way through year
      Temp = Temp[which(Data[,ReadI+1]!=-999),]   # Exclude rows with no read for this reader
      plot(x=Temp[,1],y=Temp[,2],ylim=c(0,MaxAge),xlim=c(0,MaxAge),col=rgb(red=0,green=0,blue=0,alpha=0.2),xlab="Mode predicted age | parameters",ylab="Read age",lwd=2,main=Main,pch=21,cex=0.2)
      lines(x=c(0,MaxAge),y=c(0,MaxAge), lwd=1,lty="dashed")
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI],type="l",col="red",lwd=1)
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['SD',,ReadI],type="l",col="blue",lwd=1)
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI] + 2*ErrorAndBiasArray['SD',,ReadI],type="l",col="red",lwd=1,lty="dashed")
      lines(x=ErrorAndBiasArray['True_Age',,ReadI],y=ErrorAndBiasArray['Expected_age',,ReadI] - 2*ErrorAndBiasArray['SD',,ReadI],type="l",col="red",lwd=1,lty="dashed")
    }
    mtext(side=3,outer=TRUE, text="Reads(dot), Sd(blue), expected_read(red solid line),\n and 95% CI for expected_read(red dotted line)",line=1)
  dev.off()
  
  ## AIC
  Nll = as.numeric(scan(paste(SaveFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[11])
    Df = as.numeric(scan(paste(SaveFile,"agemat.par",sep=""),comment.char="%", what="character", quiet=TRUE)[6])
    n = sum(ifelse(Data[,-1]==-999,0,1))
    Aic = 2*Nll + 2*Df
    Aicc = Aic + 2*Df*(Df+1)/(n-Df-1) 
    Bic = 2*Nll + Df*log(n)
      
  Output = list(Aic=Aic, Aicc=Aicc, Bic=Bic)
  return(Output)
}


