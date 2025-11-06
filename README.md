
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LGPclassifiers

<!-- badges: start -->
<!-- badges: end -->

This package contains all the in-house implemented classifiers for the
LGP project. For the moment it includes [Reddy Cell of Origin
(COO)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5659841/) and
originNet refined Cell of Origin (rCOO) classifiers.

## Installation

You can install the latest released version of LGPclassifiers like so:

``` r
my_repos <- c(
  "CRAN" = "https://cran.rstudio.com/",
  "BMS RSPM" = "https://pm.rdcloud.bms.com/bms-cg-biogit-bran/latest"
)
options(repos = my_repos)
install.packages("LGPclassifiers")
```

You can install the dev version:

``` r
remotes::install_git(
  "git@github.com:bms-ips/KSR-LGPclassifiers.git",
  ref = "dev",            # pin to HEAD dev
  git = "external",
  force = TRUE,
  upgrade = "always",
  build_vignettes = FALSE
)
```

## Example

Run cohort-based Reddy COO or single-sample originNet rCOO classifiers
as follows:

``` r
library(LGPclassifiers)

# Read query TPM matrix (colSums ~ 10^6)
# All three input formats ENSGs, ENSTs & Hugo names are acceptable
query.tpm <- readRDS("ensembl91-genes.salmon-tpm file from NGS360")

# Run original implementation of Reddy COO classifier
reddy <- computeCOO(query = query.tpm, useReference = F)

# Run single-sample Reddy COO classifier on reference-normalized samples (OBSOLETE)
ss.reddy <- computeCOO(query = query.tpm, useReference = T)

# Run single-sample originNet classifier on reference-normalized samples
originNet = computeOriginNet(query = query.tpm)
```
