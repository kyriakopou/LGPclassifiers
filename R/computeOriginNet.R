#' @title originNet refined COO (rCOO)
#' @description
#' Novel COO+DZ classification using originNet model
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param ... Other parameters for computeOriginNet
#' @examples
#' \dontrun{
#' # Classify query samples with origiNet
#' subType <- computeOriginNet(query = queryMatrix)
#' }
#' @export
computeOriginNet <- function(query, id2geneName = NULL, ...) {
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

  library(glmnet)

  # Check if input is a vector and convert to matrix if necessary
  if (!is.matrix(rna_norm)) {
    rna_norm <- matrix(rna_norm, ncol = 1, dimnames = list(names(rna_norm), NULL))
  }

  # coo.model <- readRDS("../originNet/build/coo.model.rds")
  coo.model <- LGPclassifiers::coo.model
  X <- t(rna_norm)[, rownames(coo.model$beta), drop = FALSE]
  pred.coo <- predict(coo.model, newx = X)[, 1]

  # dz.model <- readRDS("../originNet/build/dzsig.model.rds")
  dz.model <- LGPclassifiers::dzsig.model
  X <- t(rna_norm)[, rownames(dz.model$beta)]
  pred.dz <- predict(dz.model, newx = X)[, 1]

  df <- data.frame(ID = rownames(X), originNet_score = pred.coo, originNet = getClassFromNanostring(pred.coo),
    DZSig_score = pred.dz, DZSig = getClassFromDZSig(pred.dz)) %>%
    mutate(originNet_dz = ordered(ifelse(originNet %in% c("ABC", "Unclassified"), as.character(originNet),
      ifelse(DZSig == "DZSig+", "DZSig+", "GCB")), levels = c("GCB", "DZSig+", "Unclassified", "ABC")))

  # estimate and integrate originNet confidence in the output
  # df_conf <- getOriginNetConfidence(rna_norm, df)
  # df <- df %>%
  #   left_join(df_conf)

  return(df)

}
