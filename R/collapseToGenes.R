#' @title Collapse transcript_ids or gene_ids to gene name level
#' @description
#' Collapse transcript_ids or gene_ids to gene name level
#' @param tpm matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param geneName.map data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param featureType either "gene_id", "transcript_id" (default "gene_id")
#' @importFrom dplyr .data
#' @export
collapseToGenes <- function(tpm, geneName.map = LGPclassifiers::geneName.map) {

  tpm <- tpm[order(rownames(tpm)), , drop = FALSE]

  # if first rowname[1] contains "ENSG" or "ENST" then featureType is gene_id or transcript_id respectively
  if (grepl("^ENSG", rownames(tpm)[1])) {
    featureType <- "gene_id"
  } else if (grepl("^ENST", rownames(tpm)[1])) {
    featureType <- "transcript_id"
  } else {
    stop("Row names do not appear to be Ensembl gene or transcript IDs")
  }

  # Strip Ensembl version suffixes once (vectorized)
  ids <- sub("\\.\\d+$", "", rownames(tpm))
  # Optional cheap filter (do before mapping to cut work)
  keep <- Matrix::rowSums(tpm) > 0  # falls back to base rowSums if Matrix not loaded
  if (any(!keep)) {
    tpm <- tpm[keep, , drop = FALSE]
    ids <- ids[keep]
  }

  # Map IDs -> gene_name (vector of groups)
  # Using an index join avoids big named vectors
  gid <- geneName.map[, featureType]
  gsym <- geneName.map$gene_name

  # base match is fine once; cost is negligible vs dplyr pipeline
  idx <- match(ids, gid)

  mapped <- !is.na(idx)
  if (!any(mapped)) {
    message("Note: 0 gene/transcript IDs could be mapped to gene names")
    return(tpm[0, , drop = FALSE])
  }

  groups <- gsym[idx[mapped]]                 # gene_name per kept row

  # Aggregate rows with identical gene_name in C/Fortran via rowsum()
  # rowsum reorders groups by sort(unique(groups)) by default; preserve if desired via re-order later
  tpm.hugo <- tpm[mapped, , drop = FALSE]
  out <- rowsum(tpm.hugo, group = groups, reorder = TRUE)

  # Reporting
  n_unmapped <- sum(!mapped)
  if (n_unmapped > 0) {
    message(sprintf("Note: %d gene/transcript IDs could not be mapped to gene names", n_unmapped))
  }

  return(as.matrix(out))

}
