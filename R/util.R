#' Helper function to calculate gene mean
#' @noRd
gene.mean <- function(X, mat, bm, featureType) {
  if (featureType == "gene_id") {
    apply(mat[rownames(mat) %in% bm$gene_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, mean, na.rm = T)
  } else {
    apply(mat[rownames(mat) %in% bm$transcript_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, mean, na.rm = T)
  }
}

gene.sum <- function(X, mat, bm, featureType) {
  if (featureType == "gene_id") {
    apply(mat[rownames(mat) %in% bm$gene_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, sum, na.rm = T)
  } else {
    apply(mat[rownames(mat) %in% bm$transcript_id[bm$gene_name %in% c(X)], , drop = FALSE], 2, sum, na.rm = T)
  }
}


#' Helper function to calculate geometric mean
#' @noRd
geometric.mean <- function(x) {
  v <- na.omit(x)
  exp(mean(log(v[v > 0])))
}
