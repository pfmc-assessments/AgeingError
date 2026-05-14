# Plot output

Plots age comparisons and results from the fitted Ageing Error model

## Usage

``` r
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
```

## Arguments

- Data:

  Input data matrix

- IDataSet:

  Index of the data set used in creating the filename

- MaxAge:

  Maximum estimated age

- Report:

  Results from fitting the model

- subplot:

  Vector of which plots to create.

- Nparameters:

  Number of parameters

- LogLike:

  Negative log likelihood from fitting the model

- ReaderNames:

  Vector with names of each reader, defaults to 'Reader1', 'Reader2',
  etc. if left at the default argument of `NULL`. If you pass a vector
  of strings, the vector must be the same length as `NCOL(Data) - 1`.

- Species:

  String used at beginning of the output files

- SaveDir:

  Directory for fitted model

- verbose:

  Report messages as function runs.

- ...:

  Additional arguments passed to
  [`ageing_comparison()`](http://pfmc-assessments.github.io/AgeingError/reference/ageing_comparison.md).

## Value

Returns AIC, AICc, and BIC for fitted model.

## References

Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
Quantifying age-reading error for use in fisheries stock assessments,
with application to species in Australias southern and eastern scalefish
and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.

## Author

James T. Thorson, Ian G. Taylor
