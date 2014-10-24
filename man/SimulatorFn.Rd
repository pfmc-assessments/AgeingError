\name{SimulatorFn}
\alias{SimulatorFn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate double reading data
}
\description{
  A function to generate simulated double reading data with given properties
}
\usage{
SimulatorFn(Nreaders, M, SelexForm, ErrorParams, BiasParams, SelexParams, ReadsMat, RecCv = 0.6, RecAr1 = 0.8, Amax = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nreaders}{
  The number of ageing readers
}
  \item{M}{
  True natural mortality
}
  \item{SelexForm}{
  Form of selectivity
}
  \item{ErrorParams}{
  Error type
}
  \item{BiasParams}{
  Bias type
}
  \item{SelexParams}{
  Selectivity parameters
}
  \item{ReadsMat}{
  Matrix describing number of reads per reader combination
}
  \item{RecCv}{
  CV of recruitment
}
  \item{RecAr1}{
  First-order autoregressive coefficient for recruitment
}
  \item{Amax}{
  True maximum age
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a simulated double read matrix 

}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008. Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australias southern and eastern scalefish and shark fishery. Canadian Journal of Fisheries and Aquatic Sciences 65: 1991-2005.
}
\author{
James Thorson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#for example, please run:

example(nwfscAgeingError)

#Or consult R code
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
