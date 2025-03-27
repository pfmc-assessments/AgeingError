# AgeingError development

# AgeingError 2.1.0
* Renamed some functions to have more accurate names and match the [tidyverse style guide](https://style.tidyverse.org/functions.html)
  * `CreateData()` -> `load_data()`
  * `CreateSpecs()` -> `load_specs()`
* Added new functions to help with creating the data and specs files
  * `tally_repeats()`: Tally repeated age reading combinations in a new `count` column
  * `write_files()`: Create data and model specifications files from a data frame of age readings that can be passed to `run()` or `load_data()` and `load_specs()`
* Added automated testing for the core functions


# AgeingError 2.0.2
* Switched to the new TMB version of the package in the main branch
* Renamed github repository from nwfscAgeingError to AgeingError
* Adding vignette documenting the TMB version from Paul Burch
* Cleared out some, but not all, of the files associated with the ADMB version

# AgeingError 2.0.0

* Changed package name to AgeingError
* Added TMB code that is pre-compiled in the package
* Added `run()` as a wrapper for the code provided by @puntae
* Increased the number of examples present in the package

# nwfscAgeingError 1.3.3

* Added a `NEWS.md` file to track changes to the package.
* Included TMB agemat2.cpp source code from A. E. Punt

# nwfscAgeingError 1.3.2

This version contains updates by A. E. Punt during his work with CSIRO.

# nwfscAgeingError 1.3.1

Archive from 16th of November 2021.

# nwfscAgeingError 1.0.1

The version used for the 2017 assessment cycle.

# nwfscAgeingError 1.0.0

The version used for the 2013 assessment cycle. Though the repository was not
tagged until 2015. Melissa Haltuch noted that in 2014 the code was updated and
a summary of those updates are provided below.

 - tpl file
   - 2 new options for bias estimation.
      - (5) Splines, which requires a new argument defining the knots.
      - (6) Linear interpolation between knots, which requires a new argument
        that defines the knots.
   - Changes the calculation of the best true age because the previous way
     led to strange plots.
   - Improvement to effective sample size calculation.
 - dat file
   - Read in one data set with missing entries, where there are no double
     reads. In the 2011-based code it allowed for reading in of multiple data
     sets.
   - Increased estimation speed because of how missing values are treated.
 - R files
   - Includes code to simulate data.
   - Includes a function to format input data from a csv with readers in
     columns and double reads in the rows.
   - Plot standardized outputs from the output file. These figures could still
     use some development.
   - Ability to step through multiple runs and catalog the output using the
     new `StepwiseFn()` function.
