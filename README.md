
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
remotes::install_github(
  url = "bms-ips/KSR-LGPclassifiers",
  ref = "dev"
)
```

## Example

Run Reddy COO classifier or the internal single-sample version using a
reference dataset (ssREFERENCE) as follows:

``` r
library(LGPclassifiers)

# Read query TPM matrix (colSums ~ 10^6)
# All three input formats ENSGs, ENSTs & Hugo names are acceptable
query.tpm <- readRDS("ensembl91-genes.salmon-tpm file from NGS360")

# Run original implementation of Reddy COO classifier
reddy <- computeCOO(query = query.tpm, useReference = F)

# Run internal Reddy COO classifier on reference-normalized samples
ssREFERENCE <- computeCOO(query = query.tpm, useReference = T)
```
