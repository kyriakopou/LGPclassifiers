#' @title DLBCL90_DZ Class from DLBCL90_DZ score
#' @description
#' DLBCL90_DZ Class from DLBCL90_DZ score
#' @param score Numeric vector of DLBCL90_DZ scores
#' @export
getClassFromDLBCL90_DZ <- function(score) {
  class <- cut(score, c(-Inf, -15.6, -6.3, Inf),
    labels = c("DZSig-", "DZSig-ind", "DZSig+"))

  return(class)
}
