#' @title Normalize sample based on housekeeping genes
#' @description
#' Normalize tpm sample counts based on housekeeping gene levels
#' @param tpm.mat matrix of input samples (values must be in TPM, one row per transcript/gene id)
#' @param id2geneName data frame mapping Ensembl gene ids/transcript ids to gene names (provided as internal package data)
#' @param collapse T/F whether to collapse transcripts to gene level (default T)
#' @param featureType either "gene" or "transcript" (default "gene")
#' @param reference.mean,reference.sd reference gene means or standard deviations
#' @param houseKeeping vector of housekeeping genes to normalize to (ISY1, R3HDM1, TRIM56, UBXN4, WDR55)
#' @export
ss.normalize <- function(tpm.mat,id2geneName = NULL, collapse = TRUE, featureType = "gene",
                         houseKeeping=c("ISY1","R3HDM1","TRIM56","UBXN4","WDR55"),
                         reference.mean=NULL,reference.sd=NULL,
                         toScale = FALSE) {
  
  tmp.filtered <- tpm.mat
  
  ## Collapse transcripts to gene level
  if (collapse) {
    
    sums = round(apply(tpm.mat,2,sum),digits=0)
    if(any(sums < 1000000*0.9999 | sums > 1000000*1.0001)) {
      stop('Input tpm.mat is not a TPM matrix! Check column sums (some do not add up to 1 million)')
    }
    
    if(is.null(id2geneName)) {
      stop('A data frame mapping Ensembl gene/transcript ids to gene names is required')
    }
    
    rownames(tmp.filtered) <- gsub(x=rownames(tmp.filtered),pattern = '\\.\\d+',replacement = '')
    
    bm <- id2geneName
    tmp.filtered <- tmp.filtered+1;
    
    if (featureType=="gene") {
      bm <- bm[bm$gene_id %in% rownames(tmp.filtered),]
      tmp.filtered <- subset(x=tmp.filtered,subset=rownames(tmp.filtered) %in% bm$gene_id)
      genes <- unique(bm$gene_name);
      tmp.filtered <- do.call(rbind,lapply(genes,gene.mean,tmp.filtered,bm,featureType))
      rownames(tmp.filtered) <- genes;
    } else if (featureType=="transcript") {
      bm <- bm[bm$transcript_id %in% rownames(tmp.filtered),]
      tmp.filtered <- subset(x=tmp.filtered,subset=rownames(tmp.filtered) %in% bm$transcript_id)
      genes <- unique(bm$gene_name);
      tmp.filtered <- do.call(rbind,lapply(genes,gene.mean,tmp.filtered,bm,featureType))
      rownames(tmp.filtered) <- genes;
    } else {
      stop('Not valid featureType ("gene" or "transcript")')
    }
  }
  
  # Normalize by housekeeping genes
  TPMtmp2 <- as.matrix(tmp.filtered)
  
  # Calculate housekeeping gene normalization factor
  houseKeepingExps <- TPMtmp2[rownames(TPMtmp2) %in% houseKeeping,]
  houseKeepingExps <- apply(houseKeepingExps,2,geometric.mean)
  
  TPMtmp2 <- TPMtmp2 %*% diag(1 / houseKeepingExps)
  TPMtmp2 <- log2(TPMtmp2+1)
  colnames(TPMtmp2) <- colnames(tmp.filtered)
  # Scale data by mean and sd
  if(is.null(reference.mean) || is.null(reference.sd)) {
      TPMtmp2 <- TPMtmp2[,!is.na(colSums(TPMtmp2))]
      sDev <- apply(TPMtmp2,1,sd)
      mean <- apply(TPMtmp2,1,mean)
      output <- (TPMtmp2-mean)/sDev
      # return(list(scaled=scaled,quantile=normalize.quantiles.use.target(x = scaled,target),target=target,mean=mean,sDev=sDev,ref.transcripts=bm))
    } else {
      ref.sd <- ref.sd[match(names(ref.mean),names(ref.sd))]
      TPMtmp2 <- TPMtmp2[match(names(ref.mean),rownames(TPMtmp2)),]
      scaled <- (TPMtmp2-ref.mean)/ref.sd
      output <- scaled[,!is.na(colSums(scaled))];
      # return(list(scaled=scaled,quantile=normalize.quantiles.use.target(x=scaled,target=target),target=target,mean=ref.mean,sDev=ref.sd,ref.transcripts=ref.transcripts))
  }
  return (output)
}
