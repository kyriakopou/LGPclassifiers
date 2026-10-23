#' @title Check TPM input matrix
#' @description
#' Check that the input TPM matrix columns sum to approximately 1e6
#' @param tpm.mat matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @export
checkTPMinput <- function(tpm.mat) {
  sums <- round(apply(tpm.mat, 2, sum), digits = 0)
  # Allow some tolerance in TPM column sums (default 10%)
  tol <- 0.10
  deviating <- which(abs(sums - 1e6) > tol * 1e6)
  if (length(deviating) > 0) {
    stop(sprintf("Input tpm.mat column sums deviate > %.1f%% from 1e6 (columns: %s). Check TPM inputs.",
      tol * 100, paste(colnames(tpm.mat)[deviating], collapse = ", ")))
  } else {
    return(TRUE)
  }

}
