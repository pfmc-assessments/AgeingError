\name{RunFn}
\alias{RunFn}
\title{
  Run the model
}
\description{
  Writes DAT input file and runs the Punt et al. (2008) model programmed in ADMB
}
\usage{
RunFn(Data, SigOpt, BiasOpt, NDataSets, MinAge, MaxAge, RefAge, MinusAge, PlusAge, MaxSd, MaxExpectedAge, SaveFile, EffSampleSize = 0, Intern = TRUE, AdmbFile = NULL, JustWrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{SigOpt}{
%%     ~~Describe \code{SigOpt} here~~
}
  \item{BiasOpt}{
%%     ~~Describe \code{BiasOpt} here~~
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
  \item{MinusAge}{
%%     ~~Describe \code{MinusAge} here~~
}
  \item{PlusAge}{
%%     ~~Describe \code{PlusAge} here~~
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
  \item{AdmbFile}{
%%     ~~Describe \code{AdmbFile} here~~
}
  \item{JustWrite}{
%%     ~~Describe \code{JustWrite} here~~
}
}
\details{
}
\value{
  This does not return any values, but instead generates ADMB model output in the chosen directory
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008. Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australias southern and eastern scalefish and shark fishery. Canadian Journal of Fisheries and Aquatic Sciences 65: 1991-2005.}
\author{
James Thorson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
#for example, please run:

example(nwfscAgeingError)

#Or consult R code
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
