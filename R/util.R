#' Helper function to calculate gene mean
#' @noRd
gene.mean <- function(X,tmp.filtered,bm, featureType) {
  if (featureType == "gene") {
    apply(tmp.filtered[rownames(tmp.filtered) %in% bm$gene_id[bm$gene_name %in% c(X)],,drop = FALSE],2,mean,na.rm=T)
  } else {
    apply(tmp.filtered[rownames(tmp.filtered) %in% bm$transcript_id[bm$gene_name %in% c(X)],,drop = FALSE],2,mean,na.rm=T)
  }
}

#' Helper function to calculate geometric mean
#' @noRd
geometric.mean <- function(x) {
  v <- na.omit(x)
  exp(mean(log(v[v > 0])))
}
