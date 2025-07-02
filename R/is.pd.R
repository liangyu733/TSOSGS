#' @title is Positive Definite
#'
#' @description
#' Check whether a matrix is pd or not.
#'
#' @param X A numeric matrix.
#' @return A boolean value representing it is (if TRUE) or isn't (if FALSE) positive definite for the input matrix.
#' @noRd

is.pd <- function(X) {
  if (ncol(X) != nrow(X)) {
    warning("Non-square matrix")
    ans <- FALSE
  } else {
    ev <- eigen(X)$values
    ans <- !FALSE %in% c(abs(ev) == ev)
  }
  return(ans)
}
