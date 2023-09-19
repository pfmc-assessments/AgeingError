Description
================
nwfscAgeingError
* Wrapper for tool to estimate ageing error given double-reads of hard structures (e.g., otoliths)

  <!-- badges: start -->
  [![R build status](https://github.com/pfmc-assessments/nwfscAgeingError/workflows/R-CMD-check/badge.svg)](https://github.com/pfmc-assessments/nwfscAgeingError/actions)
  <!-- badges: end -->

Instructions
=============

First, download the executable for your operating system from the releases page:
https://github.com/pfmc-assessments/nwfscAgeingError/releases

Then install the `nwfscAgeingError` package from this github site as follows:

```r
install.packages("devtools")
devtools::install_github("pfmc-assessments/nwfscAgeingError")
# Load package
library(nwfscAgeingError)

##### Run examples
# File where the Punt et al. (2008) model (pre-compiled in ADMB) resides
SourceFile <- file.path(system.file("executables",
  package = "nwfscAgeingError"), .Platform$file.sep)
# This is where all runs will be located
dir <- getwd()
DateFile <- file.path(dir, Sys.Date())
dir.create(DateFile)

##### Generate and run with an artificial dataset
example(SimulatorFn)
utils::write.csv(AgeReads,
  file = file.path(DateFile, "Simulated_data_example.csv"))

##### Format data
example(RunFn)
##### Run the model (MAY TAKE 5-10 MINUTES)
RunFn(Data = AgeReads2, SigOpt = SigOpt, KnotAges = KnotAges,
  BiasOpt = BiasOpt,
  NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge, RefAge = 10,
  MinusAge = 1, PlusAge = 30, SaveFile = DateFile, AdmbFile = SourceFile,
  EffSampleSize = 0, Intern = FALSE, JustWrite = FALSE, CallType = "shell"
)
# Plot output
PlotOutputFn(Data = AgeReads2, MaxAge = MaxAge,
  SaveFile = DateFile, PlotType = "PDF"
)

example(StepwiseFn)
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
StepwiseFn(SearchMat = SearchMat, Data = AgeReads2,
  NDataSets = 1, MinAge = MinAge, MaxAge = MaxAge,
  RefAge = 10, MaxSd = 40, MaxExpectedAge = MaxAge+10,
  SaveFile = DateFile,
  InformationCriterion = c("AIC", "AICc", "BIC")[3]
)

```

Citing this package
=============
When using this software, please cite it as:

* Thorson, J.T., Stewart, I.J., and Punt, A.E. 2012. nwfscAgeingError: a user interface in R for the Punt et al. (2008) method for calculating ageing error and imprecision. Available from: http://github.com/pfmc-assessments/nwfscAgeingError.

and also please cite the Punt et al. (2008) paper below.

Further reading
=============
The user manual (which may not include all current features) can accessed by running the following R code:
```r
library(nwfscAgeingError)
?nwfscAgeingError
```

For more details regarding development and testing of this software please see:
* Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008. Quantifying age-reading error for use in fisheries stock assessments, with application to species in Australia’s southern and eastern scalefish and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991–2005. https://doi.org/10.1139/F08-111.


