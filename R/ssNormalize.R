# library(psych)
# library(org.Hs.eg.db)
# library(preprocessCore)

#' @title Normalize sample based on housekeeping genes
#' @description
#' Normalize sample based on housekeeping gene levels
#' @param tpm.mat matrix of input samples
#' @param id2geneName reference transcripts?
#' @param collapse T/F whether to collapse transcripts to gene level (default T)
#' @param featureType default "gene"
#' @param ref.mean,ref.sd reference gene means or standard deviations
#' @param target target
#' @param houseKeeping vector of housekeeping genes to normalize to (ISY1, R3HDM1, TRIM56, UBXN4, WDR55)
#' @param toScale T/F whether to scale (default F)
#' @export
ss.normalize <- function(tpm.mat, id2geneName,
                         collapse = TRUE,
                         featureType = "gene",
                         ref.mean = NULL, ref.sd = NULL, target = NULL,
                         houseKeeping = c("ISY1", "R3HDM1", "TRIM56", "UBXN4", "WDR55"),
                         toScale = FALSE) {
  tmp.filtered <- tpm.mat

  sums <- round(apply(tpm.mat, 2, sum), digits = 0)
  if (any(sums < 1000000 * 0.9999 | sums > 1000000 * 1.0001)) {
    stop("Input tpm.mat is not a TPM matrix! Check column sums (some do not add up to 1 million)")
  }

  rownames(tmp.filtered) <- gsub(x = rownames(tmp.filtered), pattern = "\\.\\d+", replacement = "")

  ## Collapse transcripts to gene level

  if (collapse) {
    bm <- id2geneName
    tmp.filtered <- tmp.filtered + 1
    if (featureType == "gene") {
      bm <- bm[bm$gene_id %in% rownames(tmp.filtered), ]
      genes <- unique(c(bm$gene_name, houseKeeping))
      # sorting
      genes <- genes[order(genes)]
      bm <- bm[order(bm$gene_name), ]

      # Genes with 1 occurence
      aux <- table(bm$gene_name)
      genes.1 <- names(aux[aux == 1])
      genes.2 <- names(aux[aux > 1])

      tmp.filtered.axu1 <- tmp.filtered[bm$gene_id[is.element(bm$gene_name, genes.1)], ]
      genes.1.1 <- unique(bm$gene_name[is.element(bm$gene_name, genes.1)])

      # Genes with > 1 ENSGs IDs
      tmp.filtered.aux2 <- sapply(genes.2, function(x) {
        apply(tmp.filtered[unique(bm$gene_id[bm$gene_name == x]), , drop = FALSE], 2,
          mean,
          na.rm = T
        )
      })
      genes.2.1 <- unique(bm$gene_name[is.element(bm$gene_name, genes.2)])
      tmp.filtered <- rbind(tmp.filtered.axu1, t(tmp.filtered.aux2))
      rownames(tmp.filtered) <- c(genes.1.1, genes.2.1)
    } else if (featureType == "transcript") {
      bm <- bm[bm$transcript_id %in% rownames(tmp.filtered), ]
      tmp.filtered <- subset(x = tmp.filtered, subset = rownames(tmp.filtered) %in% bm$transcript_id)
      genes <- unique(bm$gene_name)
      tmp.filtered <- do.call(rbind, lapply(genes, gene.mean, tmp.filtered, bm))
      rownames(tmp.filtered) <- genes
    } else {
      stop('Not valid featureType ("gene" or "transcript")')
    }
  }

  # Normalize by housekeeping genes
  TPMtmp2 <- as.matrix(tmp.filtered)

  # Calculate housekeeping gene normalization factor
  houseKeepingExps <- TPMtmp2[rownames(TPMtmp2) %in% houseKeeping, ]
  houseKeepingExps <- apply(houseKeepingExps, 2, geometric.mean)

  TPMtmp2 <- TPMtmp2 %*% diag(1 / houseKeepingExps)
  TPMtmp2 <- log2(TPMtmp2 + 1)
  colnames(TPMtmp2) <- colnames(tmp.filtered)
  output <- TPMtmp2
  if (toScale) {
    # Scale data by mean and sd
    if (is.null(ref.mean) || is.null(ref.sd)) {
      TPMtmp2 <- TPMtmp2[, !is.na(colSums(TPMtmp2))]
      sDev <- apply(TPMtmp2, 1, sd)
      mean <- apply(TPMtmp2, 1, mean)
      output <- (TPMtmp2 - mean) / sDev
      # return(list(scaled=scaled,quantile=normalize.quantiles.use.target(x = scaled,target),target=target,mean=mean,sDev=sDev,ref.transcripts=bm))
    } else {
      ref.sd <- ref.sd[match(names(ref.mean), names(ref.sd))]
      TPMtmp2 <- TPMtmp2[match(names(ref.mean), rownames(TPMtmp2)), ]
      scaled <- (TPMtmp2 - ref.mean) / ref.sd
      output <- scaled[, !is.na(colSums(scaled))]
      # return(list(scaled=scaled,quantile=normalize.quantiles.use.target(x=scaled,target=target),target=target,mean=ref.mean,sDev=ref.sd,ref.transcripts=ref.transcripts))
    }
  }

  return(output)
}
