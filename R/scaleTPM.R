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
  # Allow some tolerance in TPM column sums (default 5%)
  tol <- 0.10
  deviating <- which(abs(sums - 1e6) > tol * 1e6)
  if (length(deviating) > 0) {
    stop(sprintf("Input tpm.mat column sums deviate > %.1f%% from 1e6 (columns: %s). Check TPM inputs.",
      tol * 100, paste(colnames(tpm.mat)[deviating], collapse = ", ")))
  }

  # Calculate housekeeping gene normalization factor
  TPMtmp2 <- as.matrix(tpm.mat)
  houseKeepingExps <- TPMtmp2[rownames(TPMtmp2) %in% houseKeeping, ]
  houseKeepingExps <- apply(houseKeepingExps, 2, geometric.mean)

  TPMtmp2 <- TPMtmp2 %*% diag(1 / houseKeepingExps)
  TPMtmp2 <- log2(TPMtmp2 + 1)
  colnames(TPMtmp2) <- colnames(tpm.mat)

  # Maybe here an optional extra step of Quantile normalization
  # but we would have to do the same for the reference
  # if (quantile == TRUE) {
  #   TPMtmp2 <- preprocessCore::normalize.quantiles(TPMtmp2)
  # }

  # Scale data by self mean and sd
  if (is.null(ref.mean) || is.null(ref.sd)) {
    # TPMtmp2 <- TPMtmp2[, !is.na(colSums(TPMtmp2))]
    sDev <- apply(TPMtmp2, 1, sd)
    mean <- apply(TPMtmp2, 1, mean)
    output <- (TPMtmp2 - mean) / sDev
  } else {
    # use the reference mean and sd
    ref.sd <- ref.sd[match(names(ref.mean), names(ref.sd))]
    TPMtmp2 <- TPMtmp2[match(names(ref.mean), rownames(TPMtmp2)), ]
    scaled <- (TPMtmp2 - ref.mean) / ref.sd
    output <- scaled[!is.na(rowSums(scaled)), ]
  }

  return(output)
}
