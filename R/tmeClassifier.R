#' @title Apply TME26 classifier
#' @description
#' Classify samples according to TME26 classifier
#' @param exprs Input standardized expression data
#' @export
#'
tmeClassifier <- function(exprs) {
  # Get coefficients
  myModel <- LGPclassifiers::modelM3.coefficients
  modelM3.genes <- names(myModel)[-1]
  modelM3.threshold <- 3.3
  intercept <- myModel[1]

  myModel <- myModel[which(names(myModel) %in% rownames(exprs))]
  scaledData <- exprs
  scores.M3 <- apply(scaledData[toupper(names(myModel)), ], 2, function(x, y) {
    sum(x * y)
  }, myModel) + intercept
  calls.M3 <- scores.M3 > modelM3.threshold
  data_pred <- data.frame(ID = colnames(scaledData), "RNAseq_M3_Scores" = scores.M3, "RNAseq_M3_Calls" = calls.M3)

  # Determine if any TME genes are missing
  inter <- length(intersect(modelM3.genes, rownames(exprs)))
  missing <- modelM3.genes[-which(modelM3.genes %in% rownames(exprs))]
  if (inter < 25) {
    print(paste("Warning: ", 25 - inter, " of 25 TME genes missing: ", missing, sep = ""))
  }

  return(data_pred)
}
