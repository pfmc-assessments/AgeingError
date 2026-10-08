# Format an ageing-error probability matrix for FIMS

Convert a matrix with true ages in rows and observed ages in columns to
the long format expected for `ageing_error` input by FIMS.

## Usage

``` r
ageing_error_matrix_to_fims(
  probabilities,
  ages,
  fleet = NA_character_,
  timing = NA_real_
)
```

## Arguments

- probabilities:

  Numeric square matrix of probabilities. Rows represent true ages and
  columns represent observed ages.

- ages:

  Age values corresponding to the rows and columns of `probabilities`.

- fleet:

  Fleet name, or `NA` to apply to every fleet with age-composition data.

- timing:

  Model year, or `NA` to use as the default for every year.

## Value

A data frame with columns required for FIMS `ageing_error` input.
