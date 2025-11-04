#' @title Normalize TPM based on housekeeping genes and self-scale or external ref scale
#' @description
#' Normalize tpm sample counts based on housekeeping genes and self-scale or external ref scale
#' @param tpm.mat matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param ref.mean,ref.sd reference gene means or standard deviations
#' @param houseKeeping vector of housekeeping genes to normalize to (ISY1, R3HDM1, TRIM56, UBXN4, WDR55)
#' @export
scaleTPM <- function(tpm.mat, # quantile = FALSE,
                     houseKeeping = c("ISY1", "R3HDM1", "TRIM56", "UBXN4", "WDR55"),
                     ref.mean = NULL, ref.sd = NULL) {


  sums <- round(apply(tpm.mat, 2, sum), digits = 0)
  # Allow some tolerance in TPM column sums (default 10%)
  tol <- 0.10
  deviating <- which(abs(sums - 1e6) > tol * 1e6)
  if (length(deviating) > 0) {
    stop(sprintf("Input tpm.mat column sums deviate > %.1f%% from 1e6 (columns: %s). Check TPM inputs.",
      tol * 100, paste(colnames(tpm.mat)[deviating], collapse = ", ")))
  }

  TPMtmp2 <- as.matrix(tpm.mat)

  # Calculate housekeeping gene normalization factor adding eps to avoid zeros
  gm <- function(x, eps = 0.01) exp(mean(log(pmax(x, eps))))
  hk_means <- apply(TPMtmp2[rownames(TPMtmp2) %in% houseKeeping, ], 2, gm)
  TPM_hk_norm <- sweep(TPMtmp2, 2, hk_means, "/")
  TPM_hk_norm <- log2(TPM_hk_norm + 1)
  colnames(TPM_hk_norm) <- colnames(tpm.mat)


  # Scale data by self mean and sd
  if (is.null(ref.mean) || is.null(ref.sd)) {
    # TPMtmp2 <- TPMtmp2[, !is.na(colSums(TPMtmp2))]
    sDev <- apply(TPM_hk_norm, 1, sd)
    mean <- apply(TPM_hk_norm, 1, mean)
    output <- (TPM_hk_norm - mean) / sDev
  } else {
    # use the reference mean and sd
    ref.sd <- ref.sd[match(names(ref.mean), names(ref.sd))]
    TPM_hk_norm <- TPM_hk_norm[match(names(ref.mean), rownames(TPM_hk_norm)), ]
    scaled <- (TPM_hk_norm - ref.mean) / ref.sd
    output <- scaled[!is.na(rowSums(scaled)), ]
  }

  return(output)
}
