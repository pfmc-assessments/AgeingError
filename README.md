Description
================
AgeingError

* Software to estimate ageing error given double-reads of hard structures (e.g., otoliths)

  <!-- badges: start -->
  [![call-r-cmd-check](https://github.com/r4ss/r4ss/actions/workflows/call-r-cmd-check.yml/badge.svg)](https://github.com/r4ss/r4ss/actions/workflows/call-r-cmd-check.yml)
  <!-- badges: end -->

Instructions
=============

The easiest installation method is to use a pre-compiled version on R-universe via
```r
install.packages("AgeingError", repos = c("https://noaa-fisheries-integrated-toolbox.r-universe.dev", "https://cloud.r-project.org"))
```

You can also install from GitHub using the code below which will compile the package locally
(which takes longer and requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 
for Windows users).

```r
remotes::install_github("pfmc-assessments/AgeingError")
```

See [the vignette](https://pfmc-assessments.github.io/AgeingError/articles/getting_started.html) for detailed description and example use.
