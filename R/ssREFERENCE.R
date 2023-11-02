#' @title Classify sample COO using Reddy algorithm
#' @description
#' Classify sample COO as GCB-like to ABC-like. Note that this method uses a reference
#' dataset and classifies each sample relative to the reference
#' @param query matrix of normalized query gene/transcript counts for each sample
#' @param reference.mean vector of reference means (provided in package data)
#' @param reference.sd vector of reference standard deviations (provided in package data)
#' @export
ssREFERENCE <- function(query, reference.mean = ref.mean, reference.sd = ref.sd) {
  ssREFERENCE_output <- NULL

  # Get reference scaling values
  sel.genes <- intersect(names(reference.mean), rownames(query))
  reference.sd <- reference.sd[sel.genes]
  reference.mean <- reference.mean[sel.genes]

  # Scale query dataset with respect to the reference
  query.q <- (query[sel.genes, ] - reference.mean) / reference.sd

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
