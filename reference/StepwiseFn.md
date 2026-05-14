# Step-wise model selection

Run step-wise model selection to facilitate the exploration of several
modelling configurations using Akaike information criterion (AIC).

## Usage

``` r
StepwiseFn(
  SearchMat,
  Data,
  NDataSets,
  KnotAges,
  MinAge,
  MaxAge,
  RefAge,
  MaxSd,
  MaxExpectedAge,
  SaveFile,
  EffSampleSize = 0,
  Intern = TRUE,
  InformationCriterion = c("AIC", "AICc", "BIC"),
  SelectAges = TRUE
)
```

## Arguments

- SearchMat:

  A matrix explaining stepwise model selection options. One row for each
  readers error and one row for each readers bias + 2 rows, one for
  `MinusAge`, i.e., the age where the proportion at age begins to
  decrease exponentially with decreasing age, and one for `PlusAge`,
  i.e., the age where the proportion-at-age begins to decrease
  exponentially with increasing age.

  Each element of a given row is a possible value to search across for
  that reader. So, the number of columns of `SearchMat` will be the
  maximum number of options that you want to include. Think of it as
  several vectors stacked row-wise where shorter rows are filled in with
  `NA` values. If reader two only has two options that the analyst wants
  to search over the remainder of the columns should be filled with `NA`
  values for that row.

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

- NDataSets:

  This is generally `1` and other values are not implemented.

- KnotAges:

  Ages associated with each knot. This is a necessary input for
  `SigOpt = 5` or `SigOpt = 6`.

- MinAge:

  An integer, specifying the minimum possible "true" age.

- MaxAge:

  An integer, specifying the maximum possible "true" age.

- RefAge:

  An arbitrarily chosen age from which "true" age-composition
  fixed-effects are calculated as an offset. This has no effect on the
  answer but could potentially effect estimation speed.

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

- InformationCriterion:

  A string specifying the type of information criterion that should be
  used to choose the best model. The default is to use AIC, though AIC
  corrected for small sample sizes and BIC are also available.

- SelectAges:

  A logical input specifying if the boundaries should be based on
  `MinusAge` and `PlusAge`. The default is `TRUE`.

## Details

AIC seems like an appropriate method to select among possible values for
`PlusAge`, i.e., the last row of `SearchMat`, because `PlusAge`
determines the number of estimated fixed-effect hyperparameters that are
used to define the true proportion-at-age hyperdistribution. This
hyperdistribution is in turn used as a prior when integrating across a
true age associated with each otolith. This true age, which is a latent
effect, can be interpreted as a random effect with one for each
observation. So, the use of AIC to select among parameterizations of the
fixed effects defining this hyperdistribution is customary (Pinheiro and
Bates, 2009). This was tested for sablefish, where AIC lead to a true
proportion at age that was biologically plausible.

## References

Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
Quantifying age-reading error for use in fisheries stock assessments,
with application to species in Australia's southern and eastern
scalefish and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.

Pinheiro, J.C., and Bates, D. 2009. Mixed-Effects Models in S and
S-PLUS. Springer, Germany.

## See also

- [`RunFn()`](http://pfmc-assessments.github.io/AgeingError/reference/RunFn.md)
  will run a single model, where this function runs multiple models.

- [`PlotOutputFn()`](http://pfmc-assessments.github.io/AgeingError/reference/PlotOutputFn.md)
  will help summarize the output from
  [`RunFn()`](http://pfmc-assessments.github.io/AgeingError/reference/RunFn.md).

## Author

James T. Thorson

## Examples

``` r

example(RunFn)
#> 
#> RunFn> example(SimulatorFn)
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
#> 
#> RunFn> ## Not run: 
#> RunFn> ##D utils::write.csv(AgeReads,
#> RunFn> ##D   file = file.path(getwd(), "Simulated_data_example.csv")
#> RunFn> ##D )
#> RunFn> ## End(Not run)
#> RunFn> 
#> RunFn> ##### Format data
#> RunFn> Nreaders <- ncol(AgeReads)
#> 
#> RunFn> # Change NA to -999 (which the Punt software considers missing data)
#> RunFn> AgeReads <- ifelse(is.na(AgeReads), -999, AgeReads)
#> 
#> RunFn> # Potentially eliminate rows that are only read once
#> RunFn> # These rows have no information about reading error, but are potentially
#> RunFn> # informative about latent age-structure. It is unknown whether eliminating
#> RunFn> # these rows degrades estimation of error and bias, and is currently
#> RunFn> # recommended to speed up computation
#> RunFn> if (FALSE) {
#> RunFn+   KeepRow <- ifelse(
#> RunFn+     rowSums(ifelse(AgeReads == -999, 0, 1), na.rm = TRUE) <= 1,
#> RunFn+     FALSE, TRUE
#> RunFn+   )
#> RunFn+   AgeReads <- AgeReads[KeepRow, ]
#> RunFn+ }
#> 
#> RunFn> # AgeReads2 is the correctly formatted data object
#> RunFn> AgeReads2 <- rMx(c(1, AgeReads[1, ]))
#> 
#> RunFn> # Combine duplicate rows
#> RunFn> for (RowI in 2:nrow(AgeReads)) {
#> RunFn+   DupRow <- NA
#> RunFn+   for (PreviousRowJ in 1:nrow(AgeReads2)) {
#> RunFn+     if (all(
#> RunFn+       AgeReads[RowI, 1:Nreaders] == AgeReads2[PreviousRowJ, 1:Nreaders + 1]
#> RunFn+     )) {
#> RunFn+       DupRow <- PreviousRowJ
#> RunFn+     }
#> RunFn+   }
#> RunFn+   if (is.na(DupRow)) { # Add new row to AgeReads2
#> RunFn+     AgeReads2 <- rbind(AgeReads2, c(1, AgeReads[RowI, ]))
#> RunFn+   }
#> RunFn+   if (!is.na(DupRow)) { # Increment number of samples for previous duplicate
#> RunFn+     AgeReads2[DupRow, 1] <- AgeReads2[DupRow, 1] + 1
#> RunFn+   }
#> RunFn+ }
#> 
#> RunFn> ######## Determine settings for ADMB
#> RunFn> # Define minimum and maximum ages for integral across unobserved ages
#> RunFn> MinAge <- 1
#> 
#> RunFn> MaxAge <- ceiling(max(AgeReads2[, -1]) / 10) * 10
#> 
#> RunFn> BiasOpt <- c(0, -1, 0, -3)
#> 
#> RunFn> SigOpt <- c(1, -1, 6, -3)
#> 
#> RunFn> # Necessary for SigOpt option 5 or 6
#> RunFn> KnotAges <- list(NA, NA, c(1, 10, 20, MaxAge), NA)
#> 
#> RunFn> ##### Run the model (MAY TAKE 5-10 MINUTES)
#> RunFn> ## Not run: 
#> RunFn> ##D fileloc <- file.path(tempdir(), "age")
#> RunFn> ##D dir.create(fileloc, showWarnings = FALSE)
#> RunFn> ##D RunFn(
#> RunFn> ##D   Data = AgeReads2, SigOpt = SigOpt, KnotAges = KnotAges,
#> RunFn> ##D   BiasOpt = BiasOpt,
#> RunFn> ##D   NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge, RefAge = 10,
#> RunFn> ##D   MinusAge = 1, PlusAge = 30, SaveFile = fileloc,
#> RunFn> ##D   AdmbFile = file.path(system.file("executables",
#> RunFn> ##D     package = "nwfscAgeingError"
#> RunFn> ##D   ), .Platform$file.sep),
#> RunFn> ##D   EffSampleSize = 0, Intern = FALSE, JustWrite = FALSE, CallType = "shell"
#> RunFn> ##D )
#> RunFn> ## End(Not run)
#> RunFn> 
#> RunFn> 
#> RunFn> 
if (FALSE) { # \dontrun{
##### Run the model (MAY TAKE 5-10 MINUTES)
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
# Plot output
PlotOutputFn(
  Data = AgeReads2, MaxAge = MaxAge,
  SaveFile = fileloc, PlotType = "PDF"
)
} # }

##### Stepwise selection

# Parameters
MaxAge <- ceiling(max(AgeReads2) / 10) * 10
MinAge <- 1

##### Stepwise selection
StartMinusAge <- 1
StartPlusAge <- 30

# Define matrix explaining stepwise model selection options
# One row for each reader + 2 rows for
# PlusAge (age where the proportion-at-age begins to
# decrease exponentially with increasing age) and
# MinusAge (the age where the proportion-at-age begins to
# decrease exponentially with decreasing age)
# Each element of a given row is a possible value to search
# across for that reader
SearchMat <- array(NA,
  dim = c(Nreaders * 2 + 2, 7),
  dimnames = list(
    c(
      paste("Error_Reader", 1:Nreaders),
      paste("Bias_Reader", 1:Nreaders), "MinusAge", "PlusAge"
    ),
    paste("Option", 1:7)
  )
)
# Readers 1 and 3 search across options 1-3 for ERROR
SearchMat[c(1, 3), 1:3] <- rep(1, 2) %o% c(1, 2, 3)
# Reader 2 mirrors reader 1
SearchMat[2, 1] <- -1
# Reader 4 mirrors reader 3
SearchMat[4, 1] <- -3
# Reader 1 has no BIAS
SearchMat[5, 1] <- 0
# Reader 2 mirrors reader 1
SearchMat[6, 1] <- -1
# Reader 3 search across options 0-2 for BIAS
SearchMat[7, 1:3] <- c(1, 2, 0)
# Reader 4 mirrors reader 3
SearchMat[8, 1] <- -3
# MinusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
SearchMat[9, 1:7] <- c(
  StartMinusAge,
  StartMinusAge - 10,
  StartMinusAge - 4,
  StartMinusAge - 1,
  StartMinusAge + 1,
  StartMinusAge + 4,
  StartMinusAge + 10
)
SearchMat[9, 1:7] <- ifelse(SearchMat[9, 1:7] < MinAge,
  NA, SearchMat[9, 1:7]
)
# PlusAge searches with a search kernal of -10,-4,-1,+0,+1,+4,+10
SearchMat[10, 1:7] <- c(
  StartPlusAge,
  StartPlusAge - 10,
  StartPlusAge - 4,
  StartPlusAge - 1,
  StartPlusAge + 1,
  StartPlusAge + 4,
  StartPlusAge + 10
)
SearchMat[10, 1:7] <- ifelse(SearchMat[10, 1:7] > MaxAge,
  NA, SearchMat[10, 1:7]
)

# Run model selection
# This outputs a series of files
# 1. "Stepwise - Model loop X.txt" --
#   Shows the AIC/BIC/AICc value for all different combinations
#   of parameters arising from changing one parameter at a time
#   according to SearchMat during loop X
# 2. "Stepwise - Record.txt" --
#   The Xth row of IcRecord shows the record of the
#   Information Criterion for all trials in loop X,
#   while the Xth row of StateRecord shows the current selected values
#   for all parameters at the end of loop X
# 3. Standard plots for each loop
# WARNING: One run of this stepwise model building example can take
# 8+ hours, and should be run overnight
if (FALSE) { # \dontrun{
StepwiseFn(
  SearchMat = SearchMat, Data = AgeReads2,
  NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge,
  RefAge = 10, MaxSd = 40, MaxExpectedAge = MaxAge + 10,
  SaveFile = fileloc, InformationCriterion = c("AIC", "AICc", "BIC")[3]
)
} # }
```
