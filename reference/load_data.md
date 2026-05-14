# Load the formatted ageing error data

Load the formatted ageing error data

## Usage

``` r
load_data(DataFile = "data.dat", NDataSet = 1, verbose = FALSE, EchoFile = "")
```

## Arguments

- DataFile:

  Filename for input data

- NDataSet:

  Number of data sets within `DataFile`

- verbose:

  Return messages to the console (in addition to any output to
  `EchoFile`)

- EchoFile:

  A file path to a file that will be created or appended to if it
  already exists to store information about your data inputs. The
  default is `''`, which leads to output being printed to the screen
  rather than saved in a file. An example of a user-defined input would
  be `'EchoTMB.out'`.

## Author

Andre E. Punt
