#' @title COO Class from Nanostring score
#' @description
#' COO Class from Nanostring score
#' @param score Numeric vector of Nanostring scores
#' @param ... Other parameters
#' @export
getClassFromNanostring <- function(score) {
  class <- ordered(ifelse(score >= 2433.5, "ABC",
    ifelse(score <= 1907.8, "GCB", "Unclassified")), levels = c("GCB", "Unclassified", "ABC"))

  return(class)
}
