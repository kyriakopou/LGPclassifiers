#' Transcript/gene name map
#'
#' @format A data frame with 199169 rows and 3 variables:
#' \describe{
#'   \item{transcript_id}{Ensembl transcript ID}
#'   \item{gene_name}{Gene name}
#'   \item{gene_id}{Ensembl gene ID}
#' }
"geneName.map"

#' ROBUST reference distribution means
#'
#' The means of the ROBUST reference RNAseq dataset
#'
#' @format A named vector of 56478 reference gene means (after tpm scaling)
"robust.mean.tpm"

#' ROBUST reference distribution standard deviations
#'
#' The standard deviations of the ROBUST reference RNAseq dataset
#'
#' @format A named vector of 56478 reference gene sds (after tpm scaling)
"robust.sd.tpm"

#' TME26 model coefficients
#' @format A named vector of 26 coefficients, including the intercept
"modelM3.coefficients"
