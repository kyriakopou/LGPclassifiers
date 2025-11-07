#' @title Classify Sample's TME
#' @description
#' Classify matrix of RNAseq sample counts by TME26
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param ... Other parameters for scaleTPM
#' @examples
#' \dontrun{
#' # Classify query samples relative to reference dataset
#' subType <- computeTME(query = queryMatrix)
#' }
#' @export
computeTME <- function(query, id2geneName = NULL, ...) {
  # Message about using internal gene name map if not provided
  if (is.null(id2geneName)) {
    id2geneName <- LGPclassifiers::geneName.map
    message("id2geneName is NULL. Using internal gene/transcript ID to produce gene name TPMs.")
  }

  # collapse to gene-level if rownames are Ensembl gene/transcript IDs
  firstRow <- if (!is.null(rownames(query)) && length(rownames(query)) > 0) rownames(query)[1] else ""
  if (grepl("^(ENSG|ENST)", firstRow)) {
    query <- collapseToGenes(query, id2geneName)
  }


  # Use reference-based single sample classifier
  query.ss <- scaleTPM(query, ...)
  return(tmeClassifier(query.ss))
}
