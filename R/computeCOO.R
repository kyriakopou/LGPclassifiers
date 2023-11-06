#' @title Normalize sample based on housekeeping genes and run Reddy COO classifier
#' @description
#' Normalize tpm sample counts based on housekeeping gene levels
#' @param query matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param reference matrix of reference input samples (values must be in TPM, one row per transcript/gene id)
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param collapse T/F whether to collapse transcripts to gene level (default T)
#' @param featureType either "gene" or "transcript" (default "gene")
#' @param ... Other parameters for ss.normalize
#' @export
computeCOO <- function(query, reference = NULL,
                         id2geneName = NULL, collapse = TRUE, featureType = "gene", ...) {
  
  if (is.null(reference)) {
    query.ss <- ss.normalize(query,
                             id2geneName = id2geneName, collapse = collapse, featureType = featureType,
                             ...)
    return(runReddyCOO(query.ss))
  } else {
    reference.ss <- ss.normalize(reference,
                                 id2geneName = id2geneName, collapse = collapse, featureType = featureType,
                                 ...)
    ref.mean <- apply(reference.ss, 1, mean)
    ref.sd <- apply(reference.ss, 1, sd)
    query.ss <- ss.normalize(query,
                             id2geneName = id2geneName, collapse = collapse, featureType = featureType,
                             reference.mean = ref.mean, reference.sd = ref.sd,
                             ...)
    return(runReddyCOO(query.ss))
  }
}