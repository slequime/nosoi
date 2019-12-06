nosoi <img src="man/figures/logo.png" align="right" alt="" width="120" />
===============
[![Travis-CI Build Status](https://api.travis-ci.org/slequime/nosoi.svg?branch=master)](https://travis-ci.org/slequime/nosoi)
[![codecov](https://codecov.io/gh/slequime/nosoi/branch/master/graph/badge.svg)](https://codecov.io/gh/slequime/nosoi)
[![](https://img.shields.io/badge/docs-vignettes-blue.svg)](http://slequime.github.io/nosoi/)
[![](https://img.shields.io/github/license/slequime/nosoi)](http://slequime.github.io/nosoi/)

`nosoi` (pronounced no.si) is a flexible agent-based stochastic transmission chain/epidemic simulator in R, named after the *daimones* of plague, sickness and disease that escaped Pandora's jar in Greek mythology. `nosoi` is able to take into account the impact of multiple variables on the transmission process (e.g. dual-host systems such as arboviruses, within-host viral dynamics, transportation, population structure, ...), alone or taken together, to create complex but relatively intuitive epidemiological simulations.

## Installation
Stable version will be deposited on CRAN.

To get the latest (and possibly unstable) version, you can use the [`devtools`](https://github.com/hadley/devtools) package:
```R
install.packages("devtools")
devtools::install_github(repo = "slequime/nosoi")
```

## Documentation

You can find package documentation, with reference, tutorials and examples here: http://slequime.github.io/nosoi/ (built with [`pkgdown`](https://github.com/hadley/pkgdown)).
