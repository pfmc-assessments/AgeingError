

#################
#
# File structure
#
#################

# Install package
install.packages("devtools")
library("devtools")
install_github("nwfsc-assess/nwfscAgeingError")

# Load package
library(nwfscAgeingError)

# File where the Punt et al. (2008) model (pre-compiled in ADMB) resides
SourceFile = paste(system.file("executables", package="nwfscAgeingError"),"/",sep="")
  
# This is where all runs will be located
DateFile = paste(getwd(),'/',Sys.Date(),'/',sep='')
  dir.create(DateFile)

#############################
#
# Generate and run with an artificial dataset
#
# SimulatorFn() generates the data and has the following settings
#   Nreaders is the number of readers
#   ReadsMat is a matrix where each row specifies how many reads (in the first column) have a particular pattern of double reads (in the second through Nreaders+1 columns)
#   SelexForm is the selectivity-at-age form (logistic selex-at-age is the only one that is implemented)
#   SelexParams are standard to the logistic equation
#   BiasParams b in the following equation:   E[AgeRead] = b*TrueAge
#   ErrorParams are CV in the following equation:  Var[AgeRead] = (CV*TrueAge)^2
#   RecCV is the is the CV in recruitment (and recruitment is assumed stationary over time)
#   RecAR1 is first-order autoregressive coefficient in recruitment
#   Amax is the maximum allowable age
#
############################# 

##### Parameters for generating data
# This represents 2 unique readers
  # Row 1 -- Otoliths read only once by reader 
  # Row 2 -- Otoliths read twice by reader 1
  # Row 2 -- Otoliths read only once by reader 2
  # Row 4 -- Otoliths read twice by reader 2
  # Row 5 -- Otoliths read once by reader 1 and once by reader 2
Nreaders = 4
ReadsMat = cbind(NumberOfReads=rep(25,5), Reader1=c(1,1,0,0,1), Reader1_DoubleReads=c(0,1,0,0,0), Reader2=c(0,0,1,1,1), Reader2_DoubleReads=c(0,0,0,1,0))
  rownames(ReadsMat) = c("Reader1_Only", "Reader1_DoubleReads", "Reader2_Only", "Reader2_DoubleReads", "Reader1_&_Reader2")

# Generate data
AgeReads = SimulatorFn(Nreaders=Nreaders, M=0.2, SelexForm="Logistic", SelexParams=c(5,0.2), BiasParams=c(1,1,1.1,1.1), ErrorParams=c(0.2,0.2,0.2,0.2), ReadsMat=ReadsMat, RecCv=0.6, RecAr1=0.8, Amax=100)
  write.csv(AgeReads, file=paste(DateFile,"Simulated_data_example.csv",sep=""))

##### Format data  
Nreaders = ncol(AgeReads)
  AgeReads = ifelse(is.na(AgeReads),-999,AgeReads)  # Change NA to -999 (which the Punt software considers missing data)
# Potentially eliminate rows that are only read once 
  # These rows have no information about reading error, but are potentially informative about latent age-structure
  # It is unknown whether eliminating these rows degrades estimation of error and bias, and is currently recommended to speed up computation
#KeepRow = ifelse(rowSums(ifelse(AgeReads==-999,0,1),na.rm=TRUE)<=1,FALSE,TRUE)
#AgeReads = AgeReads[KeepRow,]
# Combine duplicate rows
AgeReads2 = rMx(c(1, AgeReads[1,])) # AgeReads2 is the correctly formatted data object                                  
for(RowI in 2:nrow(AgeReads)){
  DupRow = NA
  for(PreviousRowJ in 1:nrow(AgeReads2)){
    if(all(AgeReads[RowI,1:Nreaders]==AgeReads2[PreviousRowJ,1:Nreaders+1])) DupRow = PreviousRowJ
  }
  if(is.na(DupRow)) AgeReads2 = rbind(AgeReads2, c(1, AgeReads[RowI,])) # Add new row to AgeReads2
  if(!is.na(DupRow)) AgeReads2[DupRow,1] = AgeReads2[DupRow,1] + 1 # Increment number of samples for the previous duplicate
}

######## Determine settings for ADMB

# Define minimum and maximum ages for integral across unobserved ages
MinAge = 1
MaxAge = ceiling(max(AgeReads2[,-1])/10)*10

# Generate vector with settings for Bias 
# One entry for each reader
# -X = Mirror the parameters for reader X
# 0 = Unbiased (at least 1 reader has to be)
# 1 = Linear bias
# 2 = Curvilinear bias (3 param)
BiasOpt = c(0,-1,0,-3)

# Generate vector with settings for SD 
# One entry for each reader
# -X = Mirror the parameters for reader X
# 0 = No error
# 1 = Constant coefficient of variation
# 2 = Curvilinear standard deviation (3 param)
# 3 = Curvilinear coefficient of variation (3 param)
# 5 = Spline with estimated slope at beginning and end (Number of params = 2 + number of knots)
# 6 = Linear interpolation (1st knot must be 1 and last knot must be MaxAge)
SigOpt = c(1,-1,6,-3)
KnotAges = list(NA, NA, c(1,10,20,MaxAge), NA)  # Necessary for option 5 or 6

# Run the model (MAY TAKE 5-10 MINUTES)
  #Data=AgeReads2; SigOpt=SigOpt; KnotAges=KnotAges; BiasOpt=BiasOpt; NDataSets=1; MinAge=MinAge; MaxAge=MaxAge; RefAge=10; MinusAge=1; PlusAge=30; SaveFile=DateFile; AdmbFile=SourceFile; EffSampleSize=0; Intern=FALSE; JustWrite=FALSE; CallType="shell"
RunFn(Data=AgeReads2, SigOpt=SigOpt, KnotAges=KnotAges, BiasOpt=BiasOpt, NDataSets=1, MinAge=MinAge, MaxAge=MaxAge, RefAge=10, MinusAge=1, PlusAge=30, SaveFile=DateFile, AdmbFile=SourceFile, EffSampleSize=0, Intern=FALSE, JustWrite=FALSE, CallType="shell")
                 
# Plot output
  #Data=AgeReads2; MaxAge=MaxAge; SaveFile=DateFile; PlotType="PDF"
PlotOutputFn(Data=AgeReads2, MaxAge=MaxAge, SaveFile=DateFile, PlotType="PDF")

####################
#
# Stepwise selection 
#
####################

# Parameters
MaxAge = ceiling(max(AgeReads2)/10)*10
MinAge = 1

##### Stepwise selection
StartMinusAge = 1
StartPlusAge = 30

# Define matrix explaining stepwise model selection options
# One row for each reader + 2 rows for PlusAge (the age where the proportion-at-age begins to decrease exponentially with increasing age) and MinusAge (the age where the proportion-at-age begins to decrease exponentially with decreasing age)
# Each element of a given row is a possible value to search across for that reader 
SearchMat = array(NA, dim=c(Nreaders*2+2,7), dimnames=list(c(paste("Error_Reader",1:Nreaders),paste("Bias_Reader",1:Nreaders),"MinusAge","PlusAge"),paste("Option",1:7)))
  # Readers 1 and 3 search across options 1-3 for ERROR
  SearchMat[c(1,3),1:3] = rep(1,2) %o% c(1,2,3)
  # Reader 2 mirrors reader 1
  SearchMat[2,1] = -1 
  # Reader 4 mirrors reader 3
  SearchMat[4,1] = -3 
  # Reader 1 has no BIAS
  SearchMat[5,1] = 0
  # Reader 2 mirrors reader 1
  SearchMat[6,1] = -1 
  # Reader 3 search across options 0-2 for BIAS
  SearchMat[7,1:3] = c(1,2,0)
  # Reader 4 mirrors reader 3
  SearchMat[8,1] = -3 
  # MinusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
  SearchMat[9,1:7] = c(StartMinusAge,StartMinusAge-10,StartMinusAge-4,StartMinusAge-1,StartMinusAge+1,StartMinusAge+4,StartMinusAge+10)
    SearchMat[9,1:7] = ifelse(SearchMat[9,1:7]<MinAge,NA,SearchMat[9,1:7])
  # PlusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
  SearchMat[10,1:7] = c(StartPlusAge,StartPlusAge-10,StartPlusAge-4,StartPlusAge-1,StartPlusAge+1,StartPlusAge+4,StartPlusAge+10)
    SearchMat[10,1:7] = ifelse(SearchMat[10,1:7]>MaxAge,NA,SearchMat[10,1:7])
    
# Run model selection
# This outputs a series of files
  # 1. "Stepwise - Model loop X.txt" -- Shows the AIC/BIC/AICc value for all different combinations of parameters arising from changing one parameter at a time according to SearchMat during loop X
  # 2. "Stepwise - Record.txt" -- The Xth row of IcRecord shows the record of the Information Criterion for all trials in loop X, while the Xth row of StateRecord shows the current selected values for all parameters at the end of loop X
  # 3. Standard plots for each loop
# WARNING: One run of this stepwise model building example can take 8+ hours, and should be run overnight
StepwiseFn(SearchMat=SearchMat, Data=AgeReads2, NDataSets=1, MinAge=MinAge, MaxAge=MaxAge, RefAge=10, MaxSd=40, MaxExpectedAge=MaxAge+10, SaveFile=DateFile, InformationCriterion=c("AIC","AICc","BIC")[3])
  
  