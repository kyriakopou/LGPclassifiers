#' @title Classify sample COO using Reddy algorithm
#' @description
#' Classify sample COO. Note this method is implicitly Z-score based and is not "single sample".
#' It expects a "representative" distribution of COO classes for class thresholds
#' to work as expected, but scores can be used to rank from GCB-like to ABC-like regardless.
#' @param query vector of normalized query gene
#' @param ref.mean vector of reference means
#' @param ref.sd vector of reference standard deviations
#' @export
ssREFERENCE <- function(query, ref.mean = ref.mean, ref.sd = ref.sd) {
  ssREFERENCE_output <- NULL

  # Get refence scaling values
  sel.genes <- intersect(names(ref.mean), rownames(query))
  ref.sd <- ref.sd[sel.genes]
  ref.mean <- ref.mean[sel.genes]

  # Scale query dataset with respect to the reference
  query.q <- (query[sel.genes, ] - ref.mean) / ref.sd

  # Run ReddyCOO classifier
  ssREFERENCE_output <- as.data.frame(reddyCOO_new(query.q))

  return(ssREFERENCE_output)
}

#' @title Internal function calculating Reddy COO
#' @noRd
reddyCOO_new <- function(exprs) {
  # calculate reddy scores
  abc_genes <- data.frame("GeneName" = c("SH3BP5", "IRF4", "PIM1", "ENTPD1", "BLNK", "CCND2", "ETV6", "FUT8", "BMF", "IL16", "PTPN1"), "SubType" = "ABC")
  gcb_genes <- data.frame("GeneName" = c("ITPKB", "MME", "BCL6", "MYBL1", "DENND3", "NEK6", "LMO2", "LRMP", "SERPINA9"), "SubType" = "GCB")
  coo_genes <- rbind(abc_genes, gcb_genes)

  dataR_scaled <- exprs[rownames(exprs) %in% as.character(coo_genes$GeneName), ]
  # dataR_scaled <- as.data.frame(t(scale(t(dataR))))

  abc_score <- apply(subset(dataR_scaled, rownames(dataR_scaled) %in% as.character(subset(coo_genes, coo_genes$SubType == "ABC")$GeneName)), 2, mean)
  gcb_score <- apply(subset(dataR_scaled, rownames(dataR_scaled) %in% as.character(subset(coo_genes, coo_genes$SubType == "GCB")$GeneName)), 2, mean)

  final_subtype <- data.frame("RNASubtypeScore" = (abc_score - gcb_score), "Subtype" = ifelse((abc_score - gcb_score) > 0.25, "ABC", ifelse((abc_score - gcb_score) < -0.25, "GCB", "Unclassified")))
  final_subtype$Sample <- rownames(final_subtype)

  return(final_subtype)
}
