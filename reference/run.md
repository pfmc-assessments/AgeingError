# Run ageing error routine

A wrapper for running a TMB model to estimate ageing error for a given
data set and specification file.

## Usage

``` r
run(directory, file_data = "data.dat", file_specs = "data.spc")
```

## Arguments

- directory:

  A string specifying a file path to a directory where you would like to
  save the results.

- file_data:

  A string specifying the data file within 'directory'.

- file_specs:

  A string specifying the specifications file within 'directory'.

## Value

Invisibly return model output.

## See also

[`write_files()`](http://pfmc-assessments.github.io/AgeingError/reference/write_files.md)

## Author

Kelli F. Johnson

## Examples

``` r
if (FALSE) { # \dontrun{
data_test <- data.frame(
  reader1 = c(7, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, 5, 8, 5),
  reader2 = c(8, 10, 7, 6, 6, 10, 7, 9, 8, 10, 10, 5, 6, 7, 9, 7, 7, NA, NA, NA),
  reader3 = c(7, 10, 7, 6, 6, 8, 7, 9, 8, 10, 10, 5, 6, 7, NA, NA, NA, 5, 8, 5)
)
write_files(dat = data_test, dir = tempdir())
out <- run(dir = tempdir())
# see estimated parameters
out$model$par
# see model selection results
out$output$ModelSelection
# see ageing error matrices
out$output$ErrorAndBiasArray
# add to an SS3 model (assumes the model already has a single
# ageing error matrix and a maxage <= maxage in the ageing error model)
ss3_inputs <- r4ss::SS_read()
maxage <- ss3_inputs$dat$Nages
ss3_inputs$dat$ageerror <-
  out$output$ErrorAndBiasArray[c("Expected_age", "SD"), 1 + 0:maxage, "Reader 1"] |>
  as.data.frame()
r4ss::SS_write(inputlist = ss3_inputs)
} # }
```
