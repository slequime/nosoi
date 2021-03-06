nosoi <img src="man/figures/logo.png" align="right" alt="" width="120" />
===============
<!-- badges: start -->
[![Travis-CI Build Status](https://api.travis-ci.org/slequime/nosoi.svg?branch=master)](https://travis-ci.org/slequime/nosoi)
[![R build status](https://github.com/slequime/nosoi/workflows/R-CMD-check/badge.svg)](https://github.com/slequime/nosoi/actions)
[![codecov](https://codecov.io/gh/slequime/nosoi/branch/master/graph/badge.svg)](https://codecov.io/gh/slequime/nosoi)
[![](https://img.shields.io/github/license/slequime/nosoi)](http://slequime.github.io/nosoi/)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/nosoi)](https://cran.r-project.org/package=nosoi)
<!-- badges: end -->

The aim of `nosoi` (pronounced no.si) is to provide a flexible agent-based stochastic transmission chain/epidemic simulator ([Lequime et al. Methods in Ecology and Evolution 11:1002-1007](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13422)). It is named after the *daimones* of plague, sickness and disease that escaped Pandora's jar in the Greek mythology. `nosoi` is able to take into account the influence of multiple variable on the transmission process (e.g. dual-host systems (such as arboviruses), within-host viral dynamics, transportation, population structure), alone or taken together, to create complex but relatively intuitive epidemiological simulations.

## Installation
To get the current released version from CRAN:
```R
install.packages("nosoi")
```

To get the latest (and possibly unstable) version, you can use the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools")
devtools::install_github(repo = "slequime/nosoi")
```

## Documentation

You can find package documentation, with reference, tutorials and examples here: http://slequime.github.io/nosoi/ (built with [`pkgdown`](https://github.com/hadley/pkgdown)).

## Help

A Slack has been set up to help you get started or get some insight from other users. [Here is the link to join](https://join.slack.com/t/nosoi/shared_invite/zt-kk3mfk7m-d9vbb5Ze7amnR54eSvdoCw) (last updated 04/01/2021, if it does not work, just ask!).
