# Run the ageing error optimization routine

Run the ageing error optimization routine

## Usage

``` r
DoApplyAgeError(
  Species = "AgeingError",
  DataSpecs,
  ModelSpecsInp,
  AprobWght = 1e-06,
  SlopeWght = 0.01,
  SaveDir = getwd(),
  verbose = FALSE
)
```

## Arguments

- Species:

  A string that will be used to create file names. Typically, users will
  use the common name for the species of interest, especially if you are
  saving files from multiple species in a single directory. Though, the
  default is `"AgeingError"`.

- DataSpecs:

  A data object returned from
  [`load_data()`](http://pfmc-assessments.github.io/AgeingError/reference/load_data.md).

- ModelSpecsInp:

  A specification object returned from
  [`load_specs()`](http://pfmc-assessments.github.io/AgeingError/reference/load_specs.md).

- AprobWght, SlopeWght:

  Numeric values passed to the model. The defaults are 1e-06 and 0.01.
  Andre originally had these hard coded from his workspace. TODO: decide
  if they should be passed in the specifications or data files.

- SaveDir:

  A path, relative or absolute, to a directory where the results will be
  saved. The directory need not exist currently as it will be created
  dynamically.

- verbose:

  A logical specifying if messages should be printed. The default is to
  **NOT** print, i.e., `verbose = FALSE`.

## Author

Andre E. Punt
