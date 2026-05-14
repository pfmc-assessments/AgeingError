# Minimize the negative log likelihood

Minimize the negative log likelihood using `"nlmimb"` and/or `"optim"`.

## Usage

``` r
minimizer(
  model,
  method = c("optim", "nlmimb", "both"),
  lower,
  upper,
  verbose = FALSE
)
```

## Arguments

- model:

  A model to be optimized.

- method:

  A string specifying the desired method to be used for the optimization
  routine. The options are listed in the function call, where the
  default is to use `"optim"`. Using both routines is an option, via
  `"both"`, and will lead to first optimizing the model using `"nlminb"`
  and then re-optimization of the model with `"optim"`. Note that when
  using [`stats::optim()`](https://rdrr.io/r/stats/optim.html), the
  `"L-BFGS-B"` method is used rather than the default method of
  `"Nelder-Mead"`.

- lower, upper:

  Vectors of parameter bounds of the same length as the number of
  parameters in the model.

- verbose:

  A logical specifying if messages should be printed. The default is to
  **NOT** print, i.e., `verbose = FALSE`.

## Author

Andre E. Punt
