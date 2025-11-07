#' @title originNet refined COO (rCOO)
#' @description
#' Novel single-sample COO+DZ classification using internally designed originNet model
#' @param query Matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' transcript/gene id). If not specified, will default to Robust reference dataset.
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @examples
#' \dontrun{
#' # Classify query samples with origiNet
#' subType <- computeOriginNet(query = queryMatrix)
#' }
#' @importFrom glmnet predict.glmnet
#' @importFrom dplyr mutate
#' @export
computeOriginNet <- function(query, id2geneName = NULL) {
  # map to gene names if query rownames are Ensembl gene/transcript IDs
  firstRow <- if (!is.null(rownames(query)) && length(rownames(query)) > 0) rownames(query)[1] else ""
  if (grepl("^(ERCC|ENSG|ENST)", firstRow)) {
    # Message about using internal gene name map if not provided as arg
    if (is.null(id2geneName)) {
      id2geneName <- LGPclassifiers::geneName.map
      message("id2geneName is NULL. Using internal gene/transcript ID to produce gene name TPMs.")
    }

    query <- collapseToGenes(query, id2geneName)
  }

  # Check if input is a vector and convert to matrix if necessary
  if (!is.matrix(query)) {
    query <- matrix(query, ncol = 1, dimnames = list(names(query), NULL))
  }

  # TPM ss quantile normalization wrt ROBUST mean sample (SS normalisation)
  if (any(!LGPclassifiers::lgp.com.genes %in% rownames(query))) {
    # find missing genes
    missing_genes <- LGPclassifiers::lgp.com.genes[!LGPclassifiers::lgp.com.genes %in% rownames(query)]
    warning(paste("The following genes are missing from the input and will be filled with zeros:", paste(missing_genes, collapse = ", ")))
    # create a matrix with missing genes filled with zeros
    missing_matrix <- matrix(0, nrow = length(missing_genes), ncol = ncol(query),
      dimnames = list(missing_genes, colnames(query)))
    # combine the original query with the missing genes matrix
    query <- rbind(query, missing_matrix)
  }
  query.tr <- query[LGPclassifiers::lgp.com.genes, , drop = FALSE]

  ref.mean <- LGPclassifiers::robust.mean.tpm.new[LGPclassifiers::lgp.com.genes]
  query.quant <- log(quantileNormalizeToRef(query.tr, ref.mean) + 1)

  # apply coo module
  coo.model <- LGPclassifiers::originNet.coo.model
  X <- t(query.quant)[, rownames(coo.model$beta), drop = FALSE]
  pred.coo <- predict(coo.model, newx = X)[, 1]

  # apply dz module
  dz.model <- LGPclassifiers::originNet.dzsig.model
  X <- t(query.quant)[, rownames(dz.model$beta)]
  pred.dz <- predict(dz.model, newx = X)[, 1]

  df <- data.frame(ID = rownames(X), originNet_score = pred.coo, originNet = getClassFromNanostring(pred.coo),
    DZSig_score = pred.dz, DZSig = getClassFromDZSig(pred.dz)) %>%
    mutate(originNet_dz = ordered(ifelse(originNet %in% c("ABC", "Unclassified"), as.character(originNet),
      ifelse(DZSig == "DZSig+", "DZSig+", "GCB")), levels = c("GCB", "DZSig+", "Unclassified", "ABC")))

  # estimate and integrate originNet confidence in the output
  # df_conf <- getOriginNetConfidence(query.quant, df)
  # df <- df %>%
  #   left_join(df_conf)

  return(df)

}
