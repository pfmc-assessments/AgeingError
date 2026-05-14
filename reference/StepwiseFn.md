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
#> Warning: ‘RunFn’ has a help file but no examples
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
#> Error: object 'AgeReads2' not found
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
#> Error: object 'Nreaders' not found
# Readers 1 and 3 search across options 1-3 for ERROR
SearchMat[c(1, 3), 1:3] <- rep(1, 2) %o% c(1, 2, 3)
#> Error: object 'SearchMat' not found
# Reader 2 mirrors reader 1
SearchMat[2, 1] <- -1
#> Error: object 'SearchMat' not found
# Reader 4 mirrors reader 3
SearchMat[4, 1] <- -3
#> Error: object 'SearchMat' not found
# Reader 1 has no BIAS
SearchMat[5, 1] <- 0
#> Error: object 'SearchMat' not found
# Reader 2 mirrors reader 1
SearchMat[6, 1] <- -1
#> Error: object 'SearchMat' not found
# Reader 3 search across options 0-2 for BIAS
SearchMat[7, 1:3] <- c(1, 2, 0)
#> Error: object 'SearchMat' not found
# Reader 4 mirrors reader 3
SearchMat[8, 1] <- -3
#> Error: object 'SearchMat' not found
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
#> Error: object 'SearchMat' not found
SearchMat[9, 1:7] <- ifelse(SearchMat[9, 1:7] < MinAge,
  NA, SearchMat[9, 1:7]
)
#> Error: object 'SearchMat' not found
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
#> Error: object 'SearchMat' not found
SearchMat[10, 1:7] <- ifelse(SearchMat[10, 1:7] > MaxAge,
  NA, SearchMat[10, 1:7]
)
#> Error: object 'SearchMat' not found

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
