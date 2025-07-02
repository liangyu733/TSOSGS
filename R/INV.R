#' @title INVerse of matrices
#'
#' @description
#' A inverse function deal with near singular matrices by add tolParInv to the diagonals.
#'
#' @param X Squared matrix.
#' @param tolParInv tolerance Parameter of Inverse
#' @noRd

INV <- function(X, tolParInv = 10^-4) {
  if (tolParInv > 10^-1 || tolParInv < 0) {
    tolParInv <- 10^-4
    warning("'tolParInv' must range in [0, 0.1]; defaulted to 10^-4")
  }
  solve(X + diag(tolParInv, nrow(X)))
}
