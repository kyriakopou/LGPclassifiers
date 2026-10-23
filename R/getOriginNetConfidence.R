#' @title Confidence probability on originNet and DZ calls
#' @description
#' Confidence probability on originNet and DZ calls
#' @param query.quant Matrix of normalised gene expression data
#' @param df Matrix of predicted and observed calls
#' @import dplyr
#' @importFrom stats predict
#' @importFrom stats glm

getOriginNetConfidence <- function(query.quant, df) {

  df.coo <- confidenceCOO(query.quant, df)
  df.out <- confidenceDZ(query.quant, df.coo)

  return(df.out)

}
