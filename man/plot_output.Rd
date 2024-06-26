% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_output.R
\name{plot_output}
\alias{plot_output}
\title{Plot output}
\usage{
plot_output(
  Data,
  IDataSet,
  MaxAge,
  Report,
  subplot = 1:3,
  Nparameters = 0,
  LogLike = 0,
  ReaderNames = NULL,
  Species = "AgeingError",
  SaveDir = getwd(),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{Data}{Input data matrix}

\item{IDataSet}{Index of the data set used in creating the filename}

\item{MaxAge}{Maximum estimated age}

\item{Report}{Results from fitting the model}

\item{subplot}{Vector of which plots to create.}

\item{Nparameters}{Number of parameters}

\item{LogLike}{Negative log likelihood from fitting the model}

\item{ReaderNames}{Vector with names of each reader, defaults to
'Reader1', 'Reader2', etc. if left at the default argument of \code{NULL}.
If you pass a vector of strings, the vector must be the same length as
\code{NCOL(Data) - 1}.}

\item{Species}{String used at beginning of the output files}

\item{SaveDir}{Directory for fitted model}

\item{verbose}{Report messages as function runs.}

\item{...}{Additional arguments passed to \code{\link[=ageing_comparison]{ageing_comparison()}}.}
}
\value{
Returns AIC, AICc, and BIC for fitted model.
}
\description{
Plots age comparisons and results from the fitted Ageing Error model
}
\references{
Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
Quantifying age-reading error for use in fisheries stock assessments,
with application to species in Australias southern and eastern scalefish
and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.
}
\author{
James T. Thorson, Ian G. Taylor
}
