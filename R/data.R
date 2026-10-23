#' Transcript/gene name map
#'
#' @format A data frame with 199169 rows and 3 variables:
#' \describe{
#'   \item{transcript_id}{Ensembl transcript ID}
#'   \item{gene_name}{Gene name}
#'   \item{gene_id}{Ensembl gene ID}
#' }
"geneName.map"

#' ROBUST reference distribution means (used in Reddy and ssREF)
#' The means of the ROBUST reference RNAseq dataset (old TPMs after HK norm)
#' @format A named vector of 56478 reference gene means (after tpm scaling)
"robust.mean.tpm"

#' ROBUST reference distribution standard deviations (used in Reddy and ssREF)
#' The standard deviations of the ROBUST reference RNAseq dataset (old TPMs after HK norm)
#' @format A named vector of 56478 reference gene sds (after tpm scaling)
"robust.sd.tpm"

#' ROBUST reference distribution means (used on originNet)
#' The means of the ROBUST reference RNAseq dataset (new TPMs)
#' @format A named vector of 56478 reference gene means (no tpm scaling)
"robust.mean.tpm.new"

#' originNet.coo model
#' @format An object of class glmnet containing coefs for originNet.coo
"originNet.coo.model"

#' originNet.dzsig model
#' @format An object of class glmnet containing coefs for originNet.dzsig
"originNet.dzsig.model"

#' lgp.com.genes
#' @format A character vector of 14K common genes in all LGP datasets
"lgp.com.genes"

#' Robust PCA needed to compute OOD for COO calls
#' @format A list containing PCA model from ROBUST and coefficients of the logistic model to obtain class probability
"probConfIntData"

#' Robust PCA needed to compute OOD for DZSig calls
#' @format A list containing PCA model from ROBUST and coefficients of the logistic model to obtain class probability
"probConfIntData.dz"

#' DHIT signature parameters for Bayesian classification
#'
#' Pre-computed statistical parameters for the DHIT (Double Hit) signature classifier.
#' Contains mean and standard deviation values for DHITsig positive and negative distributions
#' for each of the 104 signature genes used in the Bayesian classification approach.
#'
#' @format A data frame with 104 rows (genes) and 4 columns:
#' \describe{
#'   \item{DHITsig_mean}{Mean expression value for DHITsig positive samples}
#'   \item{nonDHITsig_mean}{Mean expression value for DHITsig negative samples}
#'   \item{DHITsig_sd}{Standard deviation of expression for DHITsig positive samples}
#'   \item{nonDHITsig_sd}{Standard deviation of expression for DHITsig negative samples}
#' }
#'
#' @details
#' These parameters were derived from training data and are used by the \code{\link{runDHitSig}}
#' function to compute likelihood ratios for Bayesian classification of samples as DHIT positive
#' or negative. Row names contain the gene symbols corresponding to the DHIT signature genes.
#'
#' @seealso \code{\link{runDHitSig}}, \code{\link{dhit.model}}
"dhit_parameters"

#' DHIT signature model coefficients
#'
#' Named vector containing importance scores for each gene in the DHIT (Double Hit) signature classifier.
#' These scores weight the contribution of each gene to the final DHIT classification score.
#'
#' @format A named numeric vector with 104 elements:
#' \describe{
#'   Names are gene symbols, values are importance scores (coefficients)
#' }
#'
#' @details
#' This contains the model coefficients used by \code{\link{runDHitSig}} for efficient lookup of
#' gene importance scores during classification. Higher scores indicate stronger contribution
#' to DHIT classification. This is the main model object for the DHIT classifier.
#'
#' @seealso \code{\link{runDHitSig}}, \code{\link{dhit_parameters}}
"dhit.model"
