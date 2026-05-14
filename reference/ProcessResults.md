# Process results of the ageing error estimation

Process results of the ageing error estimation

## Usage

``` r
ProcessResults(
  Species = "AgeingError",
  SaveDir = getwd(),
  CalcEff = FALSE,
  verbose = FALSE
)
```

## Arguments

- Species:

  A string that will be used to create file names. Typically, users will
  use the common name for the species of interest, especially if you are
  saving files from multiple species in a single directory. Though, the
  default is `"AgeingError"`.

- SaveDir:

  A path, relative or absolute, to a directory where the results will be
  saved. The directory need not exist currently as it will be created
  dynamically.

- CalcEff:

  Calculate effective sample sizes (TRUE/FALSE)

- verbose:

  A logical specifying if messages should be printed. The default is to
  **NOT** print, i.e., `verbose = FALSE`.

## Author

Andre E. Punt
