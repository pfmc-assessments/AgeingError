# Tally repeated age reading combinations in a new count column

AgeingError requires a table of data with columns for each reader and
rows for each age reading combination, where a `count` column indicates
how many times the age reading combination is repeated. This function
tallies the repeated values and adds an initial `count` column with the
tally.

## Usage

``` r
tally_repeats(dat)
```

## Arguments

- dat:

  Input dataframe or tibble with columns for each reader and rows for
  each specimen

## Author

Ian G. Taylor
