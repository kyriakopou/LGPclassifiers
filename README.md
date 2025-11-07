
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LGPclassifiers

<!-- badges: start -->
<!-- badges: end -->

This package aims to include all in-house implemented classifiers for
the LGP project. For the moment it includes [Reddy Cell of Origin
(COO)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5659841/) and
originNet refined Cell of Origin (rCOO) classifiers.

## Installation

You can install the latest released version of LGPclassifiers like so:

``` r
remotes::install_github("kyriakopou/LGPclassifiers")
```

## Example

Run cohort-based Reddy COO or single-sample originNet rCOO classifiers
as follows:

``` r
library(LGPclassifiers)

# Read query TPM matrix (colSums ~ 10^6)
# All three input formats ENSGs, ENSTs & Hugo names are acceptable
query.tpm <- readRDS("ensembl91-genes.salmon-tpm file from NGS360")

# Run SOTA Reddy COO classifier
reddy <- computeCOO(query = query.tpm, useReference = F)

# Run single-sample originNet classifier on reference-normalized samples
ss.rcoo = computeOriginNet(query = query.tpm)
```
