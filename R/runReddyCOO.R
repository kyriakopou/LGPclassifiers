#' @title Classify sample COO using Reddy algorithm
#' @description
#' Classify sample COO as GCB-like to ABC-like
#' @param exprs matrix of normalized exprs counts for each sample
#' @export
runReddyCOO <- function(exprs) {
  # calculate reddy scores
  abc_genes <- data.frame("GeneName" = c("SH3BP5", "IRF4", "PIM1", "ENTPD1", "BLNK", "CCND2", "ETV6", "FUT8", "BMF", "IL16", "PTPN1"), "SubType" = "ABC")
  gcb_genes <- data.frame("GeneName" = c("ITPKB", "MME", "BCL6", "MYBL1", "DENND3", "NEK6", "LMO2", "LRMP", "SERPINA9"), "SubType" = "GCB")
  coo_genes <- rbind(abc_genes, gcb_genes)

  dataR_scaled <- exprs[rownames(exprs) %in% as.character(coo_genes$GeneName), ]

  abc_score <- apply(subset(dataR_scaled, rownames(dataR_scaled) %in% as.character(subset(coo_genes, coo_genes$SubType == "ABC")$GeneName)), 2, mean)
  gcb_score <- apply(subset(dataR_scaled, rownames(dataR_scaled) %in% as.character(subset(coo_genes, coo_genes$SubType == "GCB")$GeneName)), 2, mean)

  final_subtype <- data.frame("RNASubtypeScore" = (abc_score - gcb_score), "Subtype" = ifelse((abc_score - gcb_score) > 0.25, "ABC", ifelse((abc_score - gcb_score) < -0.25, "GCB", "Unclassified")))
  final_subtype$Sample <- rownames(final_subtype)

  return(final_subtype)
}
