#' @title single-sample quantile normalization
#' @description
#' single-sample quantile normalization to reference distribution
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param ... Other parameters for computeOriginNet
#' \dontrun{
#' # quantile normalize target counts/TPMs to reference distribution counts/TPMs
#' normalized_tpm <- quantileNormalizeToRef(target, ref)
#' }
#' @export
quantileNormalizeToRef <- function(target, ref) {

  genes <- base::intersect(rownames(target), names(ref))
  target <- target[genes, ]
  ref <- ref[genes]
  # Rank the target counts
  ranked <- apply(target, 2, rank, ties.method = "min")
  # Sort the reference distribution
  sorted_reference <- sort(ref)

  # Normalize each column in the target counts matrix
  normalized <- target
  for (i in 1:ncol(target)) {
    normalized[, i] <- sorted_reference[ranked[, i]]
  }
  return(normalized)
}
