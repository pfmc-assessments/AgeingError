\name{PlotOutputFn}
\alias{PlotOutputFn}
\title{
  Plot output
}
\description{
  Plots results from the fitted Ageing Error model
}
\usage{
  PlotOutputFn(Data, MaxAge, SaveFile, PlotType = "PDF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
Input data matrix
}
  \item{MaxAge}{
Maximum estimated age
}
  \item{SaveFile}{
Directory for fitted model 
}
  \item{PlotType}{
Type of saved plots, i.e. PDF or JPG
}
}
\details{
}
\value{
Returns AIC, AICc, and BIC for fitted model
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008. Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australias southern and eastern scalefish and shark fishery. Canadian Journal of Fisheries and Aquatic Sciences 65: 1991-2005.}
\author{
James Thorson
}
\note{
}


\seealso{
}
\examples{
#for example, please run:

example(nwfscAgeingError)

#Or consult R code

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
