#' Ageing error software
#'
#' Data input and stepwise model selection in R for the
#' Punt et al. (2008) ageing error model
#' James T. Thorson, Ian Stewart, and André E. Punt
#'
#' Function name: RunFn()
#'
#' Background:
#'
#' 	The Punt et al. (2008) model calculates the likelihood of model parameters given an observed dataset that includes age reads (henceforth 'reads') provided by multiple readers for a set of otoliths. For each reader, two sets of parameters are estimated that define the standard deviation and bias of the reads provided by that reader. Specifically, the model has parameters that approximate the expected age of each read given the true age of an otolith, and the standard deviation of a normally distributed reading error given the true age of an otolith. Each of these functional forms can be either linear or curvilinear, and each is conditioned on an unobserved 'True' age for each otolith.  This 'True' age for each otolith can be considered a random-effect, and the software computes the resulting likelihood while summing across all possible discrete values for this 'True' age for each otolith. 
#'
#'   This summation across all possible values for a 'True' age for each otolith also requires a hyperdistribution representing the 'prior' probability that an otolith is any given age; this prior is parameterized using a set of hyperparameters in addition to the parameters that govern the SD and bias for each reader.  Specifically, one hyperparameter is estimated for every age between (and including) a MinusAge and a PlusAge, which are defined exogenously for every model run.  Ages above the PlusAge or below the MinusAge have a prior Proportion-at-Age defined as a loglinear deviation from the Proportion-at-Age for the PlusAge and MinusAge.  The slope of these loglinear deviations thus constitutes an additional 1 or 2 fixed effect parameter to estimate.  The 'True' proportion-at-age is then calculated from these fixed effect and log-linear slope parameters by normalizing the resulting distribution so that it sums to one.
#'
#' Necessary Inputs:
#'
#' 	Format data: Data should be formatted with unique reading records as rows and readers/labs as columns (exampling in Table 1).  Specifically, each column corresponds to a reader, readers, lab or labs with a unique reading error and bias; the Punt (2008) model allows for approximately 15 unique columns, so the number of 'readers' must be less than this.  Additionally, an additional column inserted on the left-hand side of the data matrix indicates the number of otoliths with that unique read record; this cell is generally '1', but any instances where two or more otoliths have identical reads for all readers are combined and this cell is incremented.  Any missing entries (i.e., where a reader has not read anything for a given otolith) are indicated with a '-999' in that cell. The model can be configured such that a given column (i.e. reader) has parameter values that 'mirror' the parameter values for a reader to it's left.  This can allow estimation of a model where readers within the same lab are estimated to
#'   have the same reading error and bias.  Any instance where a particular reader (or lab) provides multiple reads for a single otolith can be dealt with by creating a 2nd column for that reader, and configuring the model so that parameters for that 2nd column mirror the parameters for the 1st column for that reader. 	Select inputs: The call-function 'FnRun()' in R writes data in the necessary format and then calls the Punt (2008) model.  This model requires several inputs, which are listed and explained below:
#'
#'   Data: This is the data set as previously formatted.  If the data has multiple rows with identical reads, this will cause an error and the 'XXX.rep' file will have a properly formatted data matrix which can be cut-pasted into a 'XXX.dat' file for use. 
#'  
#'   SigOpt: This a vector with one entry for each reader (i.e. Ncol-1 entries).  Each entry specifies the functional form of reading error as a function of true age.  Possible entries include:
#'   '-1', '-2', '-3', etc: This will make this reader mirror the estimated SD from another reader to it's left.  '-1' causes it to mirror the estimated SD for the first reader, etc.  This number has to be lower than the current entry number.
#'  
#'   '1' : Constant CV, i.e., a 1 parameter linear relationship of SD with true age.
#'  
#'   '2': Curvilinear SD, i.e., a 3 parameter Hollings-form relationship of SD with true age
#'  
#'   '3': Curvilinear with CV, i.e., a 3-parameter Hollings-form relationship of CV with true age
#'  
#'   '4': No error (but potentially bias)
#'  
#'   BiasOpt: This is a vector with one entry for each reader:
#'  
#'   '-1', '-2', '-3': See SigOpt
#'  
#'   '0': Unbiased
#'  
#'   '1': Constant CV, i.e., a 1-parameter linear relationship of bias with true age
#'  
#'   '2': Curvilinear, i.e., a 2-parameter Hollings-form relationship of bias with true age
#'  
#'   NDataSets: This is generally '1' and other values are not implemented in the current R-code. 
#'  
#'   MinAge: The minimum possible 'True' age
#'  
#'   MaxAge: The maximum possible 'True' age
#'  
#'   RefAge: An arbitrarily chosen age from which 'true' age-composition fixed-effects are calculated as an offset.  This has no effect on the answer, but could potentially effect estimation speed.
#'  
#'   MinusAge: The minimum age for which an age-specific age-composition is estimated.  Ages below this MinusAge have 'true' proportion-at-age (Pa) estimated as P_a=P_MinusAge?e^(?(MinusAge-a)), where ? is an estimated log-linear trend in the 'true' proportion-at-age.  If MinusAge = MinAge, ? is not estimated.
#'  
#'   PlusAge: Identical to MinusAge except defining the age above with age-specific age-composition is not estimated.
#'  
#'   MaxSd: An upper bound on possible values for the standard deviation of reading error
#'  
#'   MaxExpectedAge: Set to MaxAge
#'  
#'   SaveFile: Directory where 'agemat.exe' is located and where all ADMB intermediate and output files should be located.
#'  
#'   EffSampleSize: Indicating whether effective sample size should be calculated.  Missing values in the data matrix will cause this to be ineffective, in which case this should be set to '0'
#'  
#'   Intern: 'TRUE' indicates that ADMB output should be displayed in R; 'FALSE' does not.
#'
#' Stepwise model selection in R
#'
#' Function name: StepwiseFn()
#'
#' Background:
#'
#' 	Stepwise model selection allows many different model configurations to be explored: in this code, I have used AIC as the metric for comparison among model structures, although BIC or other criteria could be used.  AIC seems appropriate to select among possible PlusAge values, because this parameter determines the number of estimated fixed effect hyperparameters that are used to define the true 'Proportion-at-age' hyperdistribution.  This hyperdistribution in turn is used as a 'prior' when integrating across a 'True Age' associated with each otolith.  This 'True Age' latent effect can be interpreted as a random effect (one for each observation), so the use of AIC to select among parameterizations of the fixed effects defining this hyperdistribution is customary (Pinheiro and Bates 2009).  Additionally, the use of AIC to select the value of the PlusAge parameter appears (in preliminary analysis using Sablefish ageing error data) to lead to a 'True' proportion-at-age that is biologically plausible. 
#'
#' Necessary Inputs:
#'
#' 	Format data: Same as for a single-run
#'
#'   Select inputs: Most inputs are the same as for a single-run.  However, the 'SigOpt' 'BiasOpt' and 'PlusAge' are now specified using a matrix called 'PossibleMat', which has 2*Nreaders+2 rows and as many columns as necessary.  Row #1-#Nreaders specify the SigOpt for each reader; Next are the BiasOpt for each reader, followed by the PlusAge.  The first entry in each row specifies the starting value for that parameter in the search algorithm; a value must be specified in the first column for each parameter.  Any parameter for which the search algorithm should search across possible values has other possible values in the 2nd, 3rd, and subsequent cells in that row.  An example is given in Table 2.
#'
#' Diagnostic figures in R
#'
#' Function name: PlotOutputFn()
#'
#' Background:
#'
#' 	There are many ways to visualize the results that are provided by the Punt et al. (2008) model.  Some of these allow comparison with observed data.  However, any comparison with observed data is a little problematic, as comparisons must generally be conditioned on a 'True' age that is not observed.  In place of a 'True' age, the diagnostic plots that we present generally condition on an 'Estimated' age, which is fixed as the mode of the conditional probability-at-age for each otolith.
#'
#'   Diagnostic plots include:
#'
#'   Error and bias by reader: A panel graph where each panel shows the expected and standard deviation in age reads for that reader.  This is displayed against a scatterplot of the 'Read' and 'Estimated' ages for each otolith that was read by that reader. 
#'
#'   Proportion-at-age histogram: The estimated 'Proportion-at-age' can be plotted as a histogram, and is displayed against the 'observed' distribution of read ages.   This is useful to determine if the estimated 'proportion-at-age' is generally plausible, e.g., whether it has too many ages where the estimated proportion-at-age approaches zero (which is unlikely in a composite sample with moderate effective sample sizes).  This plot can also be used as a diagnostic to confirm that AIC has selected reasonable values for the MinusAge and PlusAge parameters. 
#'
#' Necessary Inputs:
#'
#'   The plotting function reads the 'XXX.rep' and 'XXX.par' files that are located in the directory that is provided.  It also requires specifying the MaxAge and Data, as formatted and defined earlier.
#' @references
#' Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
#' Quantifying age-reading error for use in fisheries stock assessments,
#' with application to species in Australias southern and eastern scalefish and shark fishery.
#' Canadian Journal of Fisheries and Aquatic Sciences 65: 1991-2005.
#' 
#' @examples
#' \dontrun{
#'   # File for Punt et al. (2008) model (pre-compiled in ADMB)
#'   SourceFile = paste0(
#'     system.file(package='AgeingErrorPackage'),'/executables/'
#'   )
#'  
#'   # This is where all runs will be located
#'   DateFile = paste(getwd(),'/',Sys.Date(),'/',sep='')
#'     dir.create(DateFile)
#'  
#'   #############################
#'   #
#'   # Generate and run with an artificial dataset
#'   #
#'   # SimulatorFn() generates the data and has the following settings
#'   #   Nreaders is the number of readers
#'   #   ReadsMat is a matrix where each row specifies how many reads
#'   #     (in the first column) have a particular pattern of double reads
#'   #     (in the second through Nreaders+1 columns)
#'   #   SelexForm is the selectivity-at-age form
#'   #     (logistic selex-at-age is the only one that is implemented)
#'   #   SelexParams are standard to the logistic equation
#'   #   BiasParams b in the following equation:
#'   #     E[AgeRead] = b*TrueAge
#'   #   ErrorParams are CV in the following equation:
#'   #     Var[AgeRead] = (CV*TrueAge)^2
#'   #   RecCV is the is the CV in recruitment
#'   #     (and recruitment is assumed stationary over time)
#'   #   RecAR1 is first-order autoregressive coefficient in recruitment
#'   #   Amax is the maximum allowable age
#'   #
#'   #############################
#'  
#'   ##### Parameters for generating data
#'   # This represents 2 unique readers
#'     # Row 1 -- Otoliths read only once by reader
#'     # Row 2 -- Otoliths read twice by reader 1
#'     # Row 2 -- Otoliths read only once by reader 2
#'     # Row 4 -- Otoliths read twice by reader 2
#'     # Row 5 -- Otoliths read once by reader 1 and once by reader 2
#'   Nreaders = 4
#'   ReadsMat = cbind(NumberOfReads=rep(100,5),
#'     Reader1=c(1,1,0,0,1), Reader1_DoubleReads=c(0,1,0,0,0),
#'     Reader2=c(0,0,1,1,1), Reader2_DoubleReads=c(0,0,0,1,0)
#'   )
#'   rownames(ReadsMat) = c("Reader1_Only", "Reader1_DoubleReads",
#'     "Reader2_Only", "Reader2_DoubleReads", "Reader1_&_Reader2"
#'   )
#'  
#'   # Generate data
#'   AgeReads = SimulatorFn(Nreaders = Nreaders, M = 0.2,
#'     SelexForm = "Logistic", SelexParams = c(5,0.2),
#'     BiasParams = c(1,1,1.1,1.1), ErrorParams = c(0.2,0.2,0.2,0.2),
#'     ReadsMat = ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100
#'   )
#'   utils::write.csv(AgeReads,
#'     file = paste0(DateFile,"Simulated_data_example.csv")
#'   )
#'  
#'   ##### Format data
#'   Nreaders = ncol(AgeReads)
#'   AgeReads = ifelse(is.na(AgeReads),-999,AgeReads)
#'   # Change NA to -999 (which the Punt software considers missing data)
#' 
#'   # Potentially eliminate rows that are only read once
#'     # These rows have no information about reading error,
#'     # but are potentially informative about latent age-structure
#'     # It is unknown whether eliminating these rows degrades
#'     # estimation of error and bias, and is currently recommended
#'     # to speed up computation
#'   # KeepRow = ifelse(rowSums(ifelse(AgeReads==-999,0,1),na.rm=TRUE)<=1,
#'   #   FALSE,TRUE
#'   # )
#'   # AgeReads = AgeReads[KeepRow,]
#'   # Combine duplicate rows
#'   AgeReads2 = rMx(c(1, AgeReads[1,])) # correctly formatted data object
#'   for(RowI in 2:nrow(AgeReads)){
#'     DupRow = NA
#'     for(PreviousRowJ in 1:nrow(AgeReads2)){
#'       if(all(AgeReads[RowI,1:Nreaders]==AgeReads2[PreviousRowJ,1:Nreaders+1])) DupRow = PreviousRowJ
#'     }
#'     if(is.na(DupRow)) AgeReads2 = rbind(AgeReads2, c(1, AgeReads[RowI,]))
#'     # Add new row to AgeReads2
#'     if(!is.na(DupRow)) AgeReads2[DupRow,1] = AgeReads2[DupRow,1] + 1
#'     # Increment number of samples for the previous duplicate
#'   }
#'  
#'   ######## Determine settings for ADMB
#'  
#'   # Generate vector with settings for Bias
#'   # One entry for each reader
#'   # -X = Mirror the parameters for reader X
#'   # 0 = Unbiased (at least 1 reader has to be)
#'   # 1 = Linear bias
#'   # 2 = Curvilinear bias (3 param)
#'   BiasOpt = c(0,-1,2,-3)
#'  
#'   # Generate vector with settings for SD
#'   # One entry for each reader
#'   # -X = Mirror the parameters for reader X
#'   # 0 = No error
#'   # 1 = Constant coefficient of variation
#'   # 2 = Curvilinear standard deviation (3 param)
#'   # 3 = Curvilinear coefficient of variation (3 param)
#'   SigOpt = c(3,-1,3,-3)
#'  
#'   # Define minimum and maximum ages for integral across unobserved ages
#'   MinAge = 1
#'   MaxAge = ceiling(max(AgeReads2)/10)*10
#'  
#'   # Run the model
#'   # Data=AgeReads2; SigOpt=SigOpt; BiasOpt=BiasOpt;
#'   #   NDataSets=1; MinAge=MinAge; MaxAge=MaxAge;
#'   #   RefAge=10; MinusAge=1; PlusAge=30; MaxSd=40;
#'   #   MaxExpectedAge=MaxAge+10; SaveFile=DateFile;
#'   #   AdmbFile=SourceFile; EffSampleSize=0; Intern=TRUE
#'   RunFn(Data = AgeReads2, SigOpt = SigOpt, BiasOpt = BiasOpt,
#'     NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge, RefAge = 10,
#'     MinusAge = 1, PlusAge = 30, SaveFile = DateFile, AdmbFile = SourceFile,
#'     EffSampleSize = 0, Intern = FALSE, JustWrite = FALSE
#'   )
#'  
#'   # Plot output
#'   # Data = AgeReads2; MaxAge = MaxAge; SaveFile = DateFile; PlotType = "PDF"
#'   PlotOutputFn(Data = AgeReads2, MaxAge = MaxAge, SaveFile = DateFile,
#'     PlotType = "PDF"
#'   )
#' }
#'
#' @docType package
#' @name nwfscAgeingError
NULL