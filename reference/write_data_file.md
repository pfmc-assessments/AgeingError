# Write a data file for the AgeingError package

The resulting file has the format required by the
[`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md)
function. In the future, this step could be bypassed by creating the
list output from
[`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md)
directly. Only one data set is supported (NDataSet = 1). This function
is based on a subset of the RunFn() function in the older ADMB version
of this package.

## Usage

``` r
write_data_file(
  dat,
  dir = getwd(),
  file_name = "data.dat",
  minage = 0,
  maxage = NULL,
  refage = NULL,
  minusage = NULL,
  plusage = NULL
)
```

## Arguments

- dat:

  Dataframe or tibble with columns for each reader and rows for each age
  reading combination. This could either have a count column or not. If
  not,
  [`tally_repeats()`](http://pfmc-assessments.github.io/AgeingError/reference/tally_repeats.md)
  will be called to add a count column. Order your reader/lab columns
  such that similar readers/labs are located next to one another because
  columns to the right can mirror columns to their immediate left in
  terms of parameter estimates.

- dir:

  Directory where the data file will be saved.

- file_name:

  Name of the data file.

- minage:

  An integer, specifying the minimum possible "true" age.

- maxage:

  An integer, specifying the maximum possible "true" age.

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

## Value

Invisibly returns the path to the data file
(`file.path(dir, file_name)`).

## See also

[`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md),
[`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md),
[`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md),
[`tally_repeats()`](http://pfmc-assessments.github.io/AgeingError/reference/tally_repeats.md)

## Author

Ian G. Taylor, James T. Thorson, Ian J. Stewart, Andre E. Punt

## Examples

``` r
data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)
data_file <- write_data_file(data_test, dir = tempdir(), file_name = "test.dat")
#> ℹ Input 'dat' doesn't contain a column called 'count'; adding one via tally_repeats()
#> ℹ Total observations: 20
#> ℹ Aggregated unique combinations: 12
#> ℹ Range of observed ages in the data: 5 - 10
#> ℹ Number of readers: 3
#> ℹ Max age not specified; using 15 which is the multiple of 5 which is >120% of the observed maximum
#> ℹ Minus group set to the minimum observed age 5
#> ℹ Plus group set to the maximum observed age 10
#> ℹ Reference age not specified; using 7 = floor(median(c(minusage, plusage)))
#> ℹ Writing data file to /tmp/Rtmpx7Ah1e/test.dat
```
