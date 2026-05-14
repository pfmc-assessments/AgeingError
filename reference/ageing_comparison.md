# Plot comparison of double age readings

Plot with circles proportional to how many double readings fell in each
pair of coordinates

## Usage

``` r
ageing_comparison(
  xvec,
  yvec,
  scale.pts = 2,
  col.pts = grDevices::grey(0.1, alpha = 0.5),
  col.hist = grDevices::rgb(0, 0, 0.5, alpha = 0.7),
  counts = TRUE,
  maxage = NULL,
  hist = TRUE,
  hist.frac = 0.1,
  xlab = "Age reader A",
  ylab = "Age reader B",
  title = NULL,
  png = FALSE,
  filename = "ageing_comparison.png",
  SaveFile = NULL,
  verbose = TRUE
)
```

## Arguments

- xvec:

  vector of values from reader A

- yvec:

  vector of values from reader B

- scale.pts:

  Documentation needed.

- col.pts:

  color for points

- col.hist:

  color for histograms

- counts:

  include text within each bubble showing count of values?

- maxage:

  maximum age to include in the plot (doesn't yet work well)

- hist:

  include a histogram along each axis?

- hist.frac:

  maximum value of histograms as fraction of maxage

- xlab:

  label for xvec

- ylab:

  label for yvec

- title:

  Optional title to add at top of plot

- png:

  Save plot to PNG file?

- filename:

  File name for PNG file.

- SaveFile:

  directory where plot will be saved. NULL value will make it go to
  working directory.

- verbose:

  Report messages as function runs.

## Author

Ian G. Taylor
