#' @title DZSig Class from DZSig score
#' @description
#' DZSig Class from DZSig score
#' @param score Numeric vector of DZSig scores
#' @param ... Other parameters
#' @export
getClassFromDZSig <- function(score) {

  class <- cut(score, c(-Inf, -15, 5, Inf),
    labels = c("DZSig-", NA, "DZSig+"))

  return(class)

}
