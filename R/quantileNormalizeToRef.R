#' @title single-sample quantile normalization
#' @description
#' single-sample quantile normalization to reference distribution
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param ref Numeric vector of reference distribution (one value per gene)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' \dontrun{
#' # quantile normalize query counts/TPMs to reference distribution counts/TPMs
#' query.quant <- quantileNormalizeToRef(query, ref)
#' }
#' @export
quantileNormalizeToRef <- function(query, ref) {

  genes <- base::intersect(rownames(query), names(ref))
  query <- query[genes, ]
  ref <- ref[genes]
  # Rank the query counts
  ranked <- apply(query, 2, rank, ties.method = "min")
  # Sort the reference distribution
  sorted_reference <- sort(ref)

  # Normalize each column in the query counts matrix
  normalized <- query
  for (i in 1:ncol(query)) {
    normalized[, i] <- sorted_reference[ranked[, i]]
  }
  return(normalized)
}
