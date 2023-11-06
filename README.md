
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LGPclassifiers

<!-- badges: start -->
<!-- badges: end -->

This repo contains all the in-house implemented classifiers for the LGP project. 
For the moment it includes Reddy Cell Origin classifier 
(both the original implementation and a single-sample implemented version).


## Installation

You can install the latest released version of LGPclassifiers like so:

``` r
install.packages("LGPclassifiers",
  repos = c("http://pm.rdcloud.bms.com:4242/bms-cg-biogit-bran/latest")
)
```

You can install the dev version:

``` r
remotes::install_git(url = "https://biogit.pri.bms.com/KSR/LGPclassifiers.git",
                     ref = "dev")
```

## Example

This is a basic example which shows you how generate Reddy COO
classifications for a matrix of samples by first normalizing sample tpm
counts, then calculating sample COO classifications by comparing to a
normalized reference

``` r
library(LGPclassifiers)

# Get example query matrix
# rna.counts <- readRDS("/stash/results/dev/kyriakoc/DLBCL/forManuel/rna.counts.rds")
# query <- rna.counts$NDMER

# Run single sample normalization function based on housekeeping genes
query.st <- ss.normalize(query)

# Run Reddy COO classifier on normalized samples
subType <- ssREFERENCE(query.st)
```
