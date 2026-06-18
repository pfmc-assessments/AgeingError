# Write the input files (data and specifications) for the AgeingError package

Write the input files (data and specifications) for the AgeingError
package

## Usage

``` r
write_files(
  dat,
  dir = getwd(),
  file_dat = "data.dat",
  file_specs = "data.spc",
  minage = 0,
  maxage = NULL,
  refage = NULL,
  minusage = NULL,
  plusage = NULL,
  biasopt = NULL,
  sigopt = NULL,
  knotages = NULL
)
```

## Arguments

- dat:

  Dataframe or tibble with columns for each reader and rows for each age
  reading combination. This could either have `count` as the first
  column or not. If not,
  [`tally_repeats()`](http://pfmc-assessments.github.io/AgeingError/reference/tally_repeats.md)
  will be called to add the `count` column. Order your reader/lab
  columns such that similar readers/labs are located next to one another
  because columns to the right can mirror columns to their immediate
  left in terms of parameter estimates.

- dir:

  Directory where the data file will be saved.

- minage:

  An integer, specifying the minimum possible "true" age.

- maxage:

  An integer, specifying the maximum possible "true" age. If NULL, this
  will be set to the multiple of 5 which is \>120% of the observed
  maximum age.

- refage:

  An arbitrarily chosen age from which "true" age-composition
  fixed-effects are calculated as an offset. This has no effect on the
  answer but could potentially effect estimation speed. By default this
  will be set to the maxage / 4.

- minusage:

  The minimum age for which an age-specific age-composition is
  estimated. Ages below `minusage` have "true" proportion-at-age
  (\\P\_{a}\\) estimated as \$\$P_a =
  P\_{minusage}\*exp^{(\beta\*(minusage - a))}\$\$, where beta is an
  estimated log-linear trend in the "true" proportion-at-age. If
  `minusage` = `minage`, beta is not estimated.

- plusage:

  Identical to `minusage` except defining the age above with
  age-specific age composition is not estimated.

- biasopt:

  A vector with one entry for each reader specifying the type of bias
  specific to each reader. Positive values lead to estimated parameters
  and negative values are used for shared parameters between readers.
  Parameter sharing (mirroring) is common when there is more than one
  reader in a lab working together to refine their methods such that
  they have matching techniques. If NULL is passed, the default is
  `rep(0, nreaders)` which specifies that all readers are unbiased.

  Possible entries include the following:

  `-[0-9]+`

  :   Mirror the bias of another reader, where the negative integer
      corresponds to the column of the reader that is being mirrored
      minus one, e.g., `-1` causes it to mirror reader 1. Only
      lower-numbered readers can be mirrored.

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

- sigopt:

  A vector with one entry for each reader. Each entry specifies the
  functional form of reading error as a function of true age. Positive
  values lead to estimated parameters and negative values are used for
  shared parameters between readers. If NULL is passed, the default is
  `c(1, rep(-1, nreaders - 1))` which specifies a constant CV in ageing
  error which is shared among all readers.

  Possible entries include the following:

  `-[0-9]+`

  :   Mirror the standard deviation of another reader, where the
      negative integer corresponds to the column of the reader that is
      being mirrored minus one, e.g., `-1` causes it to mirror reader 1,
      for which data is stored in the second column of `Data`. This
      number must be lower than -1 times the current position in the
      vector.

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
      of parameters is 2 + number of knots. Supported when `knotages` is
      provided.

  6

  :   Linear interpolation with a first knot of 1 and a last knot of the
      maximum age, i.e., `MaxAge`. Supported when `knotages` is
      provided.

  7

  :   A linear change in the standard deviation of random age-reading
      error, \$\_a\$, with age. This option has two parameters that need
      to be specified for each pair of independent readers in the
      specifications file.

  8

  :   A linear change in the coefficient of variation of random
      age-reading error, \$CV_a\$ , with age. This option has two
      parameters that need to be specified for each pair of independent
      readers in the specifications file.

- knotages:

  A list of knot ages for each reader. This is required when
  `sigopt = 5` or `sigopt = 6` and must have one element per reader.

## Value

Invisibly returns the path to the data file
(`file.path(dir, file_name)`).

## See also

[`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md),
[`tally_repeats()`](http://pfmc-assessments.github.io/AgeingError/reference/tally_repeats.md),
[`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md)

## Author

Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt

## Examples

``` r
data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)
write_files(dat = data_test, dir = tempdir())
#> ℹ Max age not specified; using NA which is the multiple of 5 which is >120% of the observed maximum
#> ℹ Input 'dat' doesn't contain a column called 'count'; adding one via tally_repeats()
#> ℹ Total observations: 20
#> ℹ Aggregated unique combinations: 12
#> ℹ Range of observed ages in the data: 5 - 10
#> ℹ Number of readers: 3
#> ℹ Minus group set to the minimum observed age 5
#> ℹ Plus group set to the maximum observed age 10
#> ℹ Reference age not specified; using 7 = floor(median(c(minusage, plusage)))
#> ℹ Writing data file to /tmp/RtmpAxCs7p/data.dat
#> ℹ 'biasopt' not specified; settings all readers to unbiased
#> ℹ 'sigopt' not specified; settings all readers to share a constant CV parameter
#> ℹ Writing specifications file to /tmp/RtmpAxCs7p/data.spc
run(dir = tempdir())
#> 
#> Error in determine_n_sets(file_data): number_of_sets%%1 == 0 is not TRUE
```
