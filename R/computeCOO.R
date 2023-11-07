#' @title Classify Sample Cell Of Origin
#' @description
#' Classify matrix of RNAseq sample counts by Cell Origin, either original Reddy implementation
#' or reference-based single sample classification
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param useReference Whether to use an external reference to perform a single sample classification. (default T)
#' @param reference Optional matrix of reference input samples (values must be in TPM, one row per
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param collapse T/F whether to collapse transcripts to gene level (default T)
#' @param featureType either "gene" or "transcript" (default "gene")
#' @param ... Other parameters for ss.normalize
#' @examples
#' \dontrun{
#' # Classify query samples relative to reference dataset
#' subType <- computeCOO(query = queryMatrix, useReference = T)
#' }
#' @export
computeCOO <- function(query, useReference = TRUE, reference = NULL,
                       id2geneName = NULL, collapse = TRUE, featureType = "gene", ...) {
  # Message about using internal gene name map if not provided
  if (collapse & is.null(id2geneName)) {
    id2geneName <- LGPclassifiers::geneName.map
    message("id2geneName is NULL. Using internal gene/transcript ID to gene name mapping file.")
  }

  # Use reference-based single sample classifier
  if (useReference) {
    if (is.null(reference)) {
      # message("Reference dataset is missing. Using default ROBUST reference dataset.")
      query.ss <- ss.normalize(query,
        id2geneName = id2geneName, collapse = collapse, featureType = featureType,
        ref.mean = LGPclassifiers::robust.ref.mean, ref.sd = LGPclassifiers::robust.ref.sd,
        ...
      )
      return(runReddyCOO(query.ss))
    } else {
      # Normalize and use user-provided reference dataset
      reference.ss <- ss.normalize(reference,
        id2geneName = id2geneName, collapse = collapse, featureType = featureType,
        ...
      )
      ref.mean <- apply(reference.ss, 1, mean)
      ref.sd <- apply(reference.ss, 1, sd)
      query.ss <- ss.normalize(query,
        id2geneName = id2geneName, collapse = collapse, featureType = featureType,
        ref.mean = ref.mean, ref.sd = ref.sd,
        ...
      )
      return(runReddyCOO(query.ss))
    }
  } else { # Original Reddy implementation
    query.ss <- ss.normalize(query,
      id2geneName = id2geneName, collapse = collapse, featureType = featureType,
      ...
    )
    return(runReddyCOO(query.ss))
  }
}
