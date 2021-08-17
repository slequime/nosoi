## This version 1.1.0
In this version we:
- fix typos in the documentation (following review of the manuscript) and correct an internal function (getR0).

## Test environments
* macOS-latest (Github Action, devel and release)
* ubuntu 20.04 LTS (Github Action, devel and release)
* windows-latest (Github Action, oldrel and release)
* win-builder (devel and release)

## R CMD check results

* There were no ERRORs or WARNINGs or NOTEs in our GHA and win-builder tests.

* There is one unfixed WARNING on current CRAN check result on [r-release-macos-arm64](https://www.r-project.org/nosvn/R.check/r-release-macos-arm64/nosoi-00check.html).
  This is linked with a plot using ggplot2 in the vignette.
  Package ggplot2 have the same warning on [r-release-macos-arm64](https://www.r-project.org/nosvn/R.check/r-release-macos-arm64/ggplot2-00check.html).
  We could not find a way to fix this, not having access to an ARM64 macOS machine.

## Downstream dependencies
There are currently no downstream dependencies for this package.
