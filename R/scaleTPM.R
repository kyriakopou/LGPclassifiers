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
  if (any(sums < 1000000 * 0.9999 | sums > 1000000 * 1.0001)) {
    stop("Input tpm.tpm.mat is not a TPM tpm.matrix! Check column sums (some do not add up to 1 million)")
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
    output <- scaled[, !is.na(colSums(scaled))]
  }

  return(output)
}
