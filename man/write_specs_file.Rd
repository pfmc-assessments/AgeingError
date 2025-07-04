% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_specs_file.R
\name{write_specs_file}
\alias{write_specs_file}
\title{Write a specifications file for the AgeingError package}
\usage{
write_specs_file(
  dir = getwd(),
  file_name = "data.spc",
  nreaders,
  biasopt = NULL,
  sigopt = NULL,
  knotages = NULL
)
}
\arguments{
\item{dir}{Directory where the specifications file will be saved.}

\item{file_name}{Name of the specifications file.}

\item{nreaders}{An integer, specifying the number of readers.}

\item{biasopt}{A vector with one entry for each reader specifying the
type of bias specific to each reader. Positive values lead to estimated
parameters and negative values are used for shared parameters between
readers. Parameter sharing (mirroring) is common when there
is more than one reader in a lab working together to refine their methods
such that they have matching techniques.
If NULL is passed, the default is \code{rep(0, nreaders)} which specifies that
all readers are unbiased.

Possible entries include the following:
\describe{
\item{-[0-9]+}{
Mirror the bias of another reader, where the negative integer
corresponds to the column of the reader that is being mirrored
minus one, e.g., \code{-1} causes it to mirror reader 1. Only lower-numbered
readers can be mirrored.
}
\item{0}{
Unbiased, where at least one reader has to be unbiased.
}
\item{1}{
Constant coefficient of variation, i.e., a 1-parameter linear
relationship of bias with true age.
}
\item{2}{
Curvilinear, i.e., a 2-parameter Hollings-form relationship of bias
with true age.
}
}

An example entry for the situation where you have seven readers and you
assume that the first reader is unbiased, readers 2-7 have a curvilinear
bias, reader 3 shares parameters with reader 2, reader 5 shares parameters
with reader 4, and reader 7 shares parameters with reader 6 would look
like \code{c(0, 2, -2, 2, -4, 2, -6)}.}

\item{sigopt}{A vector with one entry for each reader.
Each entry specifies the functional
form of reading error as a function of true age. Positive values lead to
estimated parameters and negative values are used for shared parameters
between readers.
If NULL is passed, the default is \code{c(1, rep(-1, nreaders - 1))} which
specifies a constant CV in ageing error which is shared among all readers.

Possible entries include the following:
\describe{
\item{-\link{0-9}+}{
Mirror the standard deviation of another reader, where the negative
integer corresponds to the column of the reader that is being
mirrored minus one, e.g., \code{-1} causes it to mirror reader 1, for
which data is stored in the second column of \code{Data}. This number must
be lower than -1 times the current position in the vector.
}
\item{0}{
No error. But, there could be potential bias.
}
\item{1}{
Constant coefficient of variation, i.e., a 1-parameter linear
relationship of the standard deviation with the true age.
}
\item{2}{
Curvilinear standard deviation, i.e., a 3-parameter Hollings-form
relationship of standard deviation with true age.
}
\item{3}{
Curvilinear coefficient of variation, i.e., a 3-parameter
Hollings-form relationship of coefficient of variation with true age.
}
\item{5}{
Spline with estimated slope at beginning and end where the number of
parameters is 2 + number of knots.
}
\item{6}{
Linear interpolation with a first knot of 1 and a last knot of the
maximum age, i.e., \code{MaxAge}.
}
}}

\item{knotages}{Ages associated with each knot. This is a necessary input
for \code{sigopt = 5} or \code{sigopt = 6}. Not implemented in this function yet.}
}
\value{
Invisibly returns the path to the specifications file (\code{file.path(dir, file_name)}).
}
\description{
The resulting file has the format required by the \code{\link[=load_specs]{load_specs()}} function.
In the future, this step could be bypassed by creating the list output
from \code{\link[=load_specs]{load_specs()}} directly.
This function is based on a subset of the RunFn() function in the older ADMB
version of this package.
}
\examples{
data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)
data_file <- write_data_file(data_test, dir = tempdir(), file_name = "test.dat")
specs_file <- write_specs_file(dir = tempdir(), nreaders = 3, file_name = "test.spc")
data <- load_data(DataFile = data_file)
specs <- load_specs(SpecsFile = specs_file, DataSpecs = data)

}
\seealso{
\code{\link[=write_files]{write_files()}}, \code{\link[=write_specs_file]{write_specs_file()}}, \code{\link[=load_specs]{load_specs()}}
}
\author{
Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt
}
