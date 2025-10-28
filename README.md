
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LGPclassifiers

<!-- badges: start -->
<!-- badges: end -->

This package contains all the in-house implemented classifiers for the
LGP project. For the moment it includes Reddy Cell Origin classifier
(both the original implementation published
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5659841/) and a
single-sample implemented version).

## Installation

You can install the latest released version of LGPclassifiers like so:

``` r
install.packages("LGPclassifiers",
  repos = c("http://pm.rdcloud.bms.com:4242/bms-cg-biogit-bran/latest")
)
```

You can install the dev version:

``` r
remotes::install_git(
  url = "https://biogit.pri.bms.com/KSR/LGPclassifiers.git",
  ref = "dev"
)
```

## Example

This is a basic example which shows you how generate Reddy COO
classifications for a matrix of samples by first normalizing sample tpm
counts, then calculating sample COO classifications by comparing to a
normalized reference

``` r
library(LGPclassifiers)

# Read example query TPM matrix (colSums ~ 10^6)
# All three input formats ENSGs, ENSTs & Hugo names are acceptable
query.tpm <- readRDS("ensembl91-genes.salmon-tpm file from NGS360")

# Run internal Reddy COO classifier on reference-normalized samples
subType <- computeCOO(query = query.tpm, useReference = T)

# Run original implementation of Reddy COO classifier
subType2 <- computeCOO(query = query.tpm, useReference = F)
```
