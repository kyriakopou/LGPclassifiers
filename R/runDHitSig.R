#' Classify sample DHit Signature Status
#'
#' @title Classify sample DHit Signature Status - Matt's implementation of Ennshi's DHitSig classifier
#' @description
#' Classify sample as DHITsig positive or negative based on Matt's implementation and thresholds of Ennshi's DHitSig classifier.
#' Uses a Bayesian approach with likelihood ratios to compute DHIT signature scores from gene expression data.
#'
#' @param exprs A numeric matrix of normalized expression counts. Rows should be genes and columns should be samples.
#'   Row names should contain gene symbols that match the DHIT signature genes.
#'
#' @return A data frame with two columns:
#'   \item{DHITscore}{Numeric DHIT signature score for each sample}
#'   \item{DHITstatus}{Logical indicating DHIT positive status (TRUE if score > 4)}
#'
#' @details
#' The function implements a Bayesian classifier that computes likelihood ratios for each gene
#' in the DHIT signature. For each sample, the algorithm:
#' \enumerate{
#'   \item For each DHIT signature gene present in the expression data
#'   \item Calculates the likelihood of the expression value under DHITsig+ and DHITsig- distributions
#'   \item Computes the log10 likelihood ratio weighted by gene importance scores
#'   \item Sums weighted log-likelihood ratios across all genes to get the final score
#'   \item Classifies samples as DHIT positive if the score exceeds 4
#' }
#'
#' The classifier uses pre-computed parameters from the \code{dhit_parameters} dataset
#' (mean and standard deviation for DHITsig+ and DHITsig- classes) and gene importance
#' scores from the \code{dhit.model} dataset.
#'
#' @references
#' Ennshi Classifier for DHIT Signature Status
#'
#' @seealso \code{\link{dhit_parameters}}, \code{\link{dhit.model}}
#'
#' @examples
#' # Create example expression data
#' genes <- c("MYC", "BCL2", "TP53")  # subset of DHIT genes
#' samples <- paste0("Sample_", 1:5)
#' exprs <- matrix(rnorm(15, mean = 0, sd = 1), nrow = 3, ncol = 5)
#' rownames(exprs) <- genes
#' colnames(exprs) <- samples
#'
#' # Run DHIT classification
#' results <- runDHitSig(exprs)
#' print(results)
#'
#' @importFrom utils data
#' @export
runDHitSig <- function(exprs) {
  # Load DHIT signature parameters and model from package data
  # Load the data objects into the local environment
  utils::data("dhit_parameters", package = "LGPclassifiers", envir = environment())
  utils::data("dhit.model", package = "LGPclassifiers", envir = environment())

  dhit_params <- dhit_parameters
  dhit_importance <- dhit.model

  # Initialize score vector
  dhitSigScore <- rep(0, ncol(exprs))
  names(dhitSigScore) <- colnames(exprs)

  # Get genes that are both in expression data and DHIT signature
  signature_genes <- intersect(rownames(exprs), rownames(dhit_params))

  # Compute DHIT scores for each sample
  for (i in 1:ncol(exprs)) {
    for (gene in signature_genes) {
      # Get expression value for this gene and sample
      value <- exprs[gene, i]

      # Get parameters for this gene
      dhit_mean <- dhit_params[gene, "DHITsig_mean"]
      non_dhit_mean <- dhit_params[gene, "nonDHITsig_mean"]
      dhit_sd <- dhit_params[gene, "DHITsig_sd"]
      non_dhit_sd <- dhit_params[gene, "nonDHITsig_sd"]

      # Calculate likelihood under DHITsig+ and DHITsig- distributions
      p1 <- dnorm(value, mean = dhit_mean, sd = dhit_sd)
      p2 <- dnorm(value, mean = non_dhit_mean, sd = non_dhit_sd)

      # Add weighted log-likelihood ratio to score
      if (p2 > 0) {  # Avoid division by zero
        log_ratio <- log(p1 / p2, base = 10)
        weight <- abs(dhit_importance[gene])
        dhitSigScore[i] <- dhitSigScore[i] + (weight * log_ratio)
      }
    }
  }

  # Create result dataframe with scores and classification
  result <- data.frame(
    DHITscore = dhitSigScore,
    # DHITstatus = dhitSigScore > 4,
    DHITstatus = getClassFromDZSig(dhitSigScore),
    stringsAsFactors = FALSE
  )

  result$ID <- rownames(result)

  return(result)
}
