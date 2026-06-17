# R-based workflow for the AgeingError package

Note: for a general introduction to the AgeingError package and its
model structure, see the [Getting Started
vignette](http://pfmc-assessments.github.io/AgeingError/articles/getting_started.md).

## Helper functions in AgeingError

The AgeingError package provides a set of helper functions to create
input files and run the ageing error model from R.

The main helpers are:

- [`write_data_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_data_file.md):
  write an input data file from an R data frame or tibble.
- [`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md):
  write an input specifications file from user inputs for bias and sigma
  options.
- [`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md):
  convenience wrapper that calls
  [`write_data_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_data_file.md)
  and
  [`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md)
  together.
- [`run()`](http://pfmc-assessments.github.io/AgeingError/reference/run.md):
  load the generated files (via
  [`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md)
  and
  [`load_specs()`](http://pfmc-assessments.github.io/AgeingError/reference/load_specs.md))
  and run the TMB estimation.

Flowchart of workflow with separate calls to
[`write_data_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_data_file.md)
and
[`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md):

Flowchart of workflow where writing the files are combined in a single
call to
[`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md):

## Writing the data file

[`write_data_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_data_file.md)
creates the data file. The main required input is an R data frame or
tibble with columns for each reader and rows for each age reading
combination (including a `count` column) or each individual fish (with
repeats), in which case the
[`tally_repeats()`](http://pfmc-assessments.github.io/AgeingError/reference/tally_repeats.md)
function will be called to add the `count` column.

Defaults will be calculated for `minage`, `maxage`, `refage`,
`minusage`, and `plusage` if not provided, but these can be overridden
by the user.

``` r

library(AgeingError)

data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)

write_data_file(
  dat = data_test,
  dir = tempdir(),
  file_name = "age_error.dat"
)
```

## Writing the specifications file

[`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md)
creates the specifications file. This file defines bias options and
sigma options for each reader. See the [Getting Started
vignette](http://pfmc-assessments.github.io/AgeingError/articles/getting_started.md)
vignette for details on the options. The options are provided as a
vector of length equal to the number of readers. When called directly,
this function also requires the number of readers (`nreaders`) to be
specified. If called via
[`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md),
the number of readers will be inferred from the data file.

``` r

write_specs_file(
  dir = tempdir(),
  file_name = "age_error.spc",
  nreaders = 4,
  biasopt = c(0, 0, 0, 0),
  sigopt = c(1, -1, 7, 8)
)
```

The options that use cubic splines for uncertainty (`sigopt = 5` and
`sigopt = 6`), the `knotages` argument must be a list of numeric knot
locations with one element per reader. Use `NA` for readers that do not
need knots.

``` r

write_specs_file(
  dir = tempdir(),
  file_name = "age_error.spc",
  nreaders = 4,
  sigopt = c(5, -1, 7, 8),
  knotages = list(c(0, 3, 5, 7, 9), NA, NA, NA)
)
```

Note that the R helper functions do not provide the user with control
over the parameter lines (low, high, or initial values). If the defaults
need to be changed, they can be modified by editing the text files
directly.

## Convenience wrapper: write_files()

[`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md)
is a convenience helper that writes both the data and specs files in a
single call.

``` r

write_files(
  dat = data_test,
  dir = tempdir(),
  file_dat = "age_error.dat",
  file_specs = "age_error.spc",
  biasopt = c(0, 0, 0),
  sigopt = c(1, -1, 2)
)
```

This is equivalent to running
[`write_data_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_data_file.md)
and
[`write_specs_file()`](http://pfmc-assessments.github.io/AgeingError/reference/write_specs_file.md)
separately and is the recommended workflow for most users.

## Running the model

[`run()`](http://pfmc-assessments.github.io/AgeingError/reference/run.md)
wraps the full process of loading the generated files and executing
AgeingError.

``` r

write_files(dat = data_test, dir = tempdir())
out <- run(dir = tempdir())

# model results
out$model$par
out$output$ModelSelection
```

Note that
[`run()`](http://pfmc-assessments.github.io/AgeingError/reference/run.md)
calls
[`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md)
and
[`load_specs()`](http://pfmc-assessments.github.io/AgeingError/reference/load_specs.md)
internally, then passes the resulting objects to
[`DoApplyAgeError()`](http://pfmc-assessments.github.io/AgeingError/reference/DoApplyAgeError.md)
and
[`ProcessResults()`](http://pfmc-assessments.github.io/AgeingError/reference/ProcessResults.md).
These functions could be called directly as needed.
