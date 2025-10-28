#' @title Classify Sample Cell Of Origin
#' @description
#' Classify matrix of RNAseq sample counts by Cell Origin, either original Reddy implementation
#' or reference-based single sample classification
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param useReference Whether to use an external reference to perform a single sample classification. (default T)
#' @param reference Optional matrix of reference input samples (values must be in TPM, one row per
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param ... Other parameters for scaleTPM
#' @examples
#' \dontrun{
#' # Classify query samples relative to reference dataset
#' subType <- computeCOO(query = queryMatrix, useReference = T)
#' }
#' @export
computeCOO <- function(query, useReference = TRUE, reference = NULL,
                       id2geneName = NULL, ...) {
  # map to gene names if query rownames are Ensembl gene/transcript IDs
  firstRow <- if (!is.null(rownames(query)) && length(rownames(query)) > 0) rownames(query)[1] else ""
  if (grepl("^(ENSG|ENST)", firstRow)) {
    # Message about using internal gene name map if not provided as arg
    if (is.null(id2geneName)) {
      id2geneName <- LGPclassifiers::geneName.map
      message("id2geneName is NULL. Using internal gene/transcript ID to produce gene name TPMs.")
    }

    query <- collapseToGenes(query, id2geneName)
  }

  # Use reference-based single sample classifier
  if (useReference) {
    if (is.null(reference)) {
      message("Scaling TPM using ROBUST as reference since no other reference was provided.")
      query.ss <- scaleTPM(query,
        ref.mean = LGPclassifiers::robust.mean.tpm, ref.sd = LGPclassifiers::robust.sd.tpm,
        ...
      )
    } else {
      # Normalize and use user-provided reference dataset
      message("Scaling TPM using user-provided reference dataset.")
      reference.ss <- scaleTPM(reference,
        ...
      )
      ref.mean <- apply(reference.ss, 1, mean)
      ref.sd <- apply(reference.ss, 1, sd)
      query.ss <- scaleTPM(query,
        ref.mean = ref.mean, ref.sd = ref.sd,
        ...
      )
    }
  } else { # Original Reddy implementation
    query.ss <- scaleTPM(query,
      ...
    )
  }
  return(runReddyCOO(query.ss))
}
