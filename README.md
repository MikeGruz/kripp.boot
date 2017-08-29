[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/kripp.boot)](https://cran.r-project.org/package=kripp.boot)

# kripp.boot

## An R Package for Performing Bootstrap Replicates of Krippendorff's alpha on Intercoder Reliability Data

*Authors: Polina Proutskova and Mike Gruszczynski*
*Maintainer: Mike Gruszczynski*

This function implements Prof. Klaus Krippendorff's algorithm for bootstrapping the Krippendorff's alpha coefficient. It computes confidence values (reliability estimates) for the given probabilities. 

Because Krippendorff's alpha has no *a priori* known distribution, deriving confidence intervals from alpha estimates requires bootstrap replicates to be generated. This package allows those replicates to be generated quickly and easily.

## Installation

To install the package, type the following:

```
install.packages("kripp.boot")
library(kripp.boot)
```


Or you can install the development version from GitHub:

```
library(devtools)
install_github("mikegruz/kripp.boot")
library(kripp.boot)
```

## Usage

```
kripp.boot(x, iter = 2000, probs = c(.025, .975), 
           method = c("nominal", "ordinal", "interval", "ratio"))
```

Where:

`x` is a matrix corresponding to judges and columns corresponding to rated objects.

`iter` is the number of bootstrap replicates to perform.

`probs` is a vector of probabilities for which confidence values are computed.

`method` is the metric used to calculate the difference function in computing alpha. Methods for nominal, ordinal, interval, and ratio data are currently implemented.
