\name{StepwiseFn}
\alias{StepwiseFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Stepwise model selection
}
\description{
  Code for running stepwise model selection
}
\usage{
StepwiseFn(SearchMat, Data, NDataSets, MinAge, MaxAge, RefAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, Intern = TRUE, InformationCriterion = "AIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SearchMat}{
%%     ~~Describe \code{SearchMat} here~~
}
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{NDataSets}{
%%     ~~Describe \code{NDataSets} here~~
}
  \item{MinAge}{
%%     ~~Describe \code{MinAge} here~~
}
  \item{MaxAge}{
%%     ~~Describe \code{MaxAge} here~~
}
  \item{RefAge}{
%%     ~~Describe \code{RefAge} here~~
}
  \item{MaxSd}{
%%     ~~Describe \code{MaxSd} here~~
}
  \item{MaxExpectedAge}{
%%     ~~Describe \code{MaxExpectedAge} here~~
}
  \item{SaveFile}{
%%     ~~Describe \code{SaveFile} here~~
}
  \item{EffSampleSize}{
%%     ~~Describe \code{EffSampleSize} here~~
}
  \item{Intern}{
%%     ~~Describe \code{Intern} here~~
}
  \item{InformationCriterion}{
%%     ~~Describe \code{InformationCriterion} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008. Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australias southern and eastern scalefish and shark fishery. Canadian Journal of Fisheries and Aquatic Sciences 65: 1991-2005.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

if(FALSE){
  # Run initial model 
  example(nwfscAgeingError)
  
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
}  
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
