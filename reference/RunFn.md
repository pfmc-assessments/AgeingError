# Run ageing error model

Run the Punt et al. (2008) ADMB-based ageing error model from within R.

## Usage

``` r
RunFn(
  Data,
  SigOpt,
  KnotAges,
  BiasOpt,
  NDataSets = 1,
  MinAge,
  MaxAge,
  RefAge,
  MinusAge,
  PlusAge,
  MaxSd,
  MaxExpectedAge,
  SaveFile,
  EffSampleSize = 0,
  Intern = TRUE,
  AdmbFile = NULL,
  JustWrite = FALSE,
  CallType = "system",
  ExtraArgs = " -est",
  verbose = TRUE
)
```

## Arguments

- Data:

  This is the data set with the first column being an integer providing
  the number of otoliths that are included in the row and the subsequent
  columns are the reader or lab estimated ag,e where each reader/lab has
  a unique reading error and bias. The modeling framework allows for, at
  most, 15 readers, i.e., 16 columns. There should not be any identical
  rows in the data frame because otoliths that have the exact same read
  from every reader/lab should be combined into a single row with the
  count as the first column. If you failed to combine identical rows
  prior to running the model, you will be alerted with an error and the
  `XXX.rep` file will have a properly formatted data which can be'
  cut-pasted into a `XXX.dat` file for use. Missing reads from a given
  reader/lab should be entered as `-999`. Order your reader/lab columns
  such that similar readers/labs are located next to one another because
  columns to the right can mirror columns to their immediate left in
  terms of parameter estimates.

- SigOpt:

  This a vector with one entry for each reader (i.e.,
  `length(SigOpt) == NCOL(Data) -1`). Each entry specifies the
  functional form of reading error as a function of true age. Possible
  entries include the following:

  -0-9+

  :   Mirror the standard deviation of another reader, where the
      negative integer corresponds to the column of the reader/lab that
      is being mirrored minus one, e.g., `-1` causes it to mirror
      reader/lab 1, for which data is stored in the second column of
      `Data`. This number must be lower than -1 times the current
      position in the vector.

  0

  :   No error. But, there could be potential bias.

  1

  :   Constant coefficient of variation, i.e., a 1-parameter linear
      relationship of the standard deviation with the true age.

  2

  :   Curvilinear standard deviation, i.e., a 3-parameter Hollings-form
      relationship of standard deviation with true age.

  3

  :   Curvilinear coefficient of variation, i.e., a 3-parameter
      Hollings-form relationship of coefficient of variation with true
      age.

  5

  :   Spline with estimated slope at beginning and end where the number
      of parameters is 2 + number of knots.

  6

  :   Linear interpolation with a first knot of 1 and a last knot of the
      maximum age, i.e., `MaxAge`.

- KnotAges:

  Ages associated with each knot. This is a necessary input for
  `SigOpt = 5` or `SigOpt = 6`.

- BiasOpt:

  A vector with one entry for each reader/lab specifying the type of
  bias specific to each reader. Positive values lead to estimated
  parameters and negative values are used for shared parameters between
  readers, just like with `SigOpt`. Parameter sharing is common when
  there is more than one reader in a lab working together to refine
  their methods such that they have matching techniques. Possible
  entries include the following:

  -0-9+

  :   Mirror the bias of another reader, where the negative integer
      corresponds to the column of the reader/lab that is being mirrored
      minus one, e.g., `-1` causes it to mirror reader/lab 1, for which
      data is stored in the second column of `Data`. This number must be
      lower than -1 times the current position in the vector.

  0

  :   Unbiased, where at least one reader has to be unbiased.

  1

  :   Constant coefficient of variation, i.e., a 1-parameter linear
      relationship of bias with true age.

  2

  :   Curvilinear, i.e., a 2-parameter Hollings-form relationship of
      bias with true age.

  An example entry for the situation where you have seven readers and
  you assume that the first reader is unbiased, readers 2-7 have a
  curvilinear bias, reader 3 shares parameters with reader 2, reader 5
  shares parameters with reader 4, and reader 7 shares parameters with
  reader 6 would look like `c(0, 2, -2, 2, -4, 2, -6)`.

- NDataSets:

  This is generally `1` and other values are not implemented.

- MinAge:

  An integer, specifying the minimum possible "true" age.

- MaxAge:

  An integer, specifying the maximum possible "true" age.

- RefAge:

  An arbitrarily chosen age from which "true" age-composition
  fixed-effects are calculated as an offset. This has no effect on the
  answer but could potentially effect estimation speed.

- MinusAge:

  The minimum age for which an age-specific age-composition is
  estimated. Ages below `MinusAge` have "true" proportion-at-age
  (\\P\_{a}\\) estimated as \$\$P_a =
  P\_{MinusAge}\*exp^{(\beta\*(MinusAge - a))}\$\$, where beta is an
  estimated log-linear trend in the "true" proportion-at-age. If
  `MinusAge` = `MinAge`, beta is not estimated.

- PlusAge:

  Identical to `MinusAge` except defining the age above with
  age-specific age composition is not estimated.

- MaxSd:

  An upper bound on possible values for the standard deviation of
  reading error.

- MaxExpectedAge:

  Set to MaxAge.

- SaveFile:

  Directory where `agemat.exe` is located and where all ADMB
  intermediate and output files should be located. If `AdmbFile` is
  specified then `agemat.exe` is copied from that directory to
  `SaveFile`.

- EffSampleSize:

  Indicating whether effective sample size should be calculated. Missing
  values in the data matrix will cause this to be ineffective, in which
  case this should be set to `0`.

- Intern:

  A logical input that controls the amount of output displayed, where
  `TRUE` indicates that ADMB output should be displayed in R and `FALSE`
  leads to the suppression of this information.

- AdmbFile:

  An optional character entry that specifies the directory from which
  `agemat.exe` is to be copied from to `SaveFile`.

- JustWrite:

  A logical input that allows just the data files to be written without
  running ADMB executable.

- CallType:

  Either `"system"` or `"shell"` depending on Operating System or how R
  is being run. The default is `"system"`.

- ExtraArgs:

  A string of characters providing extra arguments passed to ADMB. The
  default is `" -est"`.

- verbose:

  A logical input that controls the amount of feedback users receive
  from the program. The default is to provide the most output as
  possible with `verbose = TRUE`.

## Details

The premise of Punt *et al.* (2008) is to calculate the likelihood of
model parameters given an observed data set of otolith age reads from
multiple age readers. For each reader/lab, two parameters are defined,
one for standard deviation and one for bias. The model calculates the
expected age of each read and the standard deviation of a normally
distributed reading error given the true age of an otolith. These
relationships can be linear or curvilinear.

The true age is obviously an unobserved process and can be considered a
random effect. Thus, the software computes the likelihood while summing
across all possible discrete values for the true age of each otolith.
This true age requires a hyperdistribution that represents the prior
probability that an otolith is any given age. The hyperdistribution is
controlled by a set of hyperparameters and the parameters that govern
the standard deviation and bias of each age reader/lab. Specifically,
one hyperparameter is estimated for every age between and including the
`MinusAge` and `PlusAge`. Ages outside of this range have a prior
proportion at age defined as a loglinear deviation from the proportion
at age for the extreme ages, i.e., `MinusAge` and `PlusAge`. The slope
of these loglinear deviations thus constitutes an additional 1 or 2
fixed effect parameters. The true proportion at age is then calculated
from these fixed effects and loglinear slope parameters by normalizing
the resulting distribution such that it sums to one.

## See also

- [`StepwiseFn()`](http://pfmc-assessments.github.io/AgeingError/reference/StepwiseFn.md)
  will run multiple models.

- [`PlotOutputFn()`](http://pfmc-assessments.github.io/AgeingError/reference/PlotOutputFn.md)
  will help summarize the output.

## Author

James T. Thorson, Ian J. Stewart, Andre E. Punt, Ian G. Taylor

## Examples

``` r
example(SimulatorFn)
#> 
#> SmltrF> # Parameters for generating data
#> SmltrF> # This represents 2 unique readers
#> SmltrF> # Row 1 -- Otoliths read only once by reader
#> SmltrF> # Row 2 -- Otoliths read twice by reader 1
#> SmltrF> # Row 2 -- Otoliths read only once by reader 2
#> SmltrF> # Row 4 -- Otoliths read twice by reader 2
#> SmltrF> # Row 5 -- Otoliths read once by reader 1 and once by reader 2
#> SmltrF> ReadsMat <- structure(matrix(
#> SmltrF+   nrow = 5, ncol = 5,
#> SmltrF+   c(
#> SmltrF+     rep(25, 5),
#> SmltrF+     1, 1, 0, 0, 1,
#> SmltrF+     0, 1, 0, 0, 0,
#> SmltrF+     0, 0, 1, 1, 1,
#> SmltrF+     0, 0, 0, 1, 0
#> SmltrF+   )
#> SmltrF+ ), dimnames = list(
#> SmltrF+   c(
#> SmltrF+     "Reader1_Only", "Reader1_DoubleReads",
#> SmltrF+     "Reader2_Only", "Reader2_DoubleReads",
#> SmltrF+     "Reader1_&_Reader2"
#> SmltrF+   ),
#> SmltrF+   c(
#> SmltrF+     "NumberOfReads",
#> SmltrF+     "Reader1", "Reader1_DoubleReads",
#> SmltrF+     "Reader2", "Reader2_DoubleReads"
#> SmltrF+   )
#> SmltrF+ ))
#> 
#> SmltrF> # Generate data
#> SmltrF> set.seed(2)
#> 
#> SmltrF> AgeReads <- SimulatorFn(
#> SmltrF+   Nreaders = 4, M = 0.2,
#> SmltrF+   SelexForm = "Logistic",
#> SmltrF+   SelexParams = c(5, 0.2), BiasParams = c(1, 1, 1.1, 1.1),
#> SmltrF+   ErrorParams = c(0.2, 0.2, 0.2, 0.2), ReadsMat = ReadsMat,
#> SmltrF+   RecCv = 0.6, RecAr1 = 0.8, Amax = 100
#> SmltrF+ )
if (FALSE) { # \dontrun{
utils::write.csv(AgeReads,
  file = file.path(getwd(), "Simulated_data_example.csv")
)
} # }

##### Format data
Nreaders <- ncol(AgeReads)
# Change NA to -999 (which the Punt software considers missing data)
AgeReads <- ifelse(is.na(AgeReads), -999, AgeReads)

# Potentially eliminate rows that are only read once
# These rows have no information about reading error, but are potentially
# informative about latent age-structure. It is unknown whether eliminating
# these rows degrades estimation of error and bias, and is currently
# recommended to speed up computation
if (FALSE) {
  KeepRow <- ifelse(
    rowSums(ifelse(AgeReads == -999, 0, 1), na.rm = TRUE) <= 1,
    FALSE, TRUE
  )
  AgeReads <- AgeReads[KeepRow, ]
}

# AgeReads2 is the correctly formatted data object
AgeReads2 <- rMx(c(1, AgeReads[1, ]))

# Combine duplicate rows
for (RowI in 2:nrow(AgeReads)) {
  DupRow <- NA
  for (PreviousRowJ in 1:nrow(AgeReads2)) {
    if (all(
      AgeReads[RowI, 1:Nreaders] == AgeReads2[PreviousRowJ, 1:Nreaders + 1]
    )) {
      DupRow <- PreviousRowJ
    }
  }
  if (is.na(DupRow)) { # Add new row to AgeReads2
    AgeReads2 <- rbind(AgeReads2, c(1, AgeReads[RowI, ]))
  }
  if (!is.na(DupRow)) { # Increment number of samples for previous duplicate
    AgeReads2[DupRow, 1] <- AgeReads2[DupRow, 1] + 1
  }
}

######## Determine settings for ADMB
# Define minimum and maximum ages for integral across unobserved ages
MinAge <- 1
MaxAge <- ceiling(max(AgeReads2[, -1]) / 10) * 10
BiasOpt <- c(0, -1, 0, -3)
SigOpt <- c(1, -1, 6, -3)
# Necessary for SigOpt option 5 or 6
KnotAges <- list(NA, NA, c(1, 10, 20, MaxAge), NA)

##### Run the model (MAY TAKE 5-10 MINUTES)
if (FALSE) { # \dontrun{
fileloc <- file.path(tempdir(), "age")
dir.create(fileloc, showWarnings = FALSE)
RunFn(
  Data = AgeReads2, SigOpt = SigOpt, KnotAges = KnotAges,
  BiasOpt = BiasOpt,
  NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge, RefAge = 10,
  MinusAge = 1, PlusAge = 30, SaveFile = fileloc,
  AdmbFile = file.path(system.file("executables",
    package = "nwfscAgeingError"
  ), .Platform$file.sep),
  EffSampleSize = 0, Intern = FALSE, JustWrite = FALSE, CallType = "shell"
)
} # }
```
