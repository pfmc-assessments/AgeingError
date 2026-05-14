# Plot output

Plots age comparisons and results from the fitted model. Comparisons
must be conditioned on a true age that is not observed. And, in place of
a true age, the diagnostic plots generally condition on an estimated
age, which is fixed as the mode of the conditional probability at age
for each otolith.

## Usage

``` r
PlotOutputFn(
  Data,
  MaxAge,
  SaveFile,
  PlotType = c("PNG", "PDF"),
  subplot = 1:3,
  ReaderNames = NULL,
  ...
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

- MaxAge:

  An integer, specifying the maximum possible "true" age.

- SaveFile:

  Directory where `agemat.exe` is located and where all ADMB
  intermediate and output files should be located. If `AdmbFile` is
  specified then `agemat.exe` is copied from that directory to
  `SaveFile`.

- PlotType:

  A string specifying the type of saved plots that you desire. The
  default is to save `.png` files via an argument of `"PNG"`. The other
  option is to save `.pdf` files via `"PDF"`.

- subplot:

  Vector of integers specifying which plots to create. The default is to
  create three plots.

- ReaderNames:

  Vector with names of each reader, defaults to "Reader 1", "Reader 2",
  etc.

- ...:

  Additional arguments passed to
  [`ageing_comparison()`](http://pfmc-assessments.github.io/AgeingError/reference/ageing_comparison.md).

## Value

Returns AIC, AICc, and BIC for fitted model.

## Details

1.  Error and bias by reader/lab: A panel graph is provided where each
    panel shows the expected and standard deviation in age reads for
    that reader/lab. This is displayed against a scatter plot of the
    read and estimated ages for each otolith that was read by that
    reader/lab.

2.  Proportion-at-age histogram: The estimated proportion at age can be
    plotted as a histogram and is displayed against the observed
    distribution of read ages. This is useful to determine if hte
    estimated proportion at age is generally plausible, e.g., whether it
    has too many ages where the estimated proportion at age approaches
    zero, which is unlikely in a composite sample with moderate
    effective sample sizes. This plot can also be used as a diagnostic
    to confirm that AIC has selected reasonable values for the
    `MinusAge` and `PlusAge` parameters.

The function will read in `XXX.rep` and `XXX.par` files that are located
in `SaveFile`.

## References

Punt, A.E., Smith, D.C., KrusicGolub, K., and Robertson, S. 2008.
Quantifying age-reading error for use in fisheries stock assessments,
with application to species in Australias southern and eastern scalefish
and shark fishery. Can. J. Fish. Aquat. Sci. 65: 1991-2005.

## See also

- [`RunFn()`](http://pfmc-assessments.github.io/AgeingError/reference/RunFn.md)

- [`StepwiseFn()`](http://pfmc-assessments.github.io/AgeingError/reference/StepwiseFn.md)

## Author

James T. Thorson, Ian G. Taylor
