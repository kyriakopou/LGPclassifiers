#' @title Classify Sample's TME
#' @description
#' Classify matrix of RNAseq sample counts by TME26
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param featureType either "gene_id" or "transcript_id" or "gene.name" (no collapse function is executed) (default "gene_id")
#' @param ... Other parameters for scaleTPM
#' @examples
#' \dontrun{
#' # Classify query samples relative to reference dataset
#' subType <- computeTME(query = queryMatrix)
#' }
#' @export
computeTME <- function(query, id2geneName = NULL, featureType = "gene_id", ...) {
  # Message about using internal gene name map if not provided
  if (is.null(id2geneName)) {
    id2geneName <- LGPclassifiers::geneName.map
    message("id2geneName is NULL. Using internal gene/transcript ID to produce gene name TPMs.")
  }

  ## Collapse transcripts to gene level if gene_ids or transcript_ids are given
  if (featureType %in% c("gene_id", "transcript_id")) {
    query <- collapseToGenes(query, id2geneName, featureType = featureType)
  }

  # Use reference-based single sample classifier
  query.ss <- scaleTPM(query, ...)
  return(tmeClassifier(query.ss))
}
