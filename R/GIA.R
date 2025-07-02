#' @title Geometrically Isomorphic Allocation for rectangular field
#'
#' @description
#' Deal with the geometrically isomorphic allocation by reflecting and rotating to standard order.
#' A 2 by 3 example is demonstrated.
#' The allocation matrix(c(1:6), 2, 3, byrow = TRUE) is equivalent to:
#' 1. matrix(c(3,2,1,6,5,4), 2, 3, byrow = TRUE) by reflection,
#' 2. matrix(c(4,5,6,1,2,3), 2, 3, byrow = TRUE) by reflection,
#' 3. matrix(c(6,5,4,3,2,1), 2, 3, byrow = TRUE) by rotation.
#' NULL if consider those are different.
#'
#' @param v is a vector
#' @param gia.rf A integer vector of length 2 indicating nrow and ncol of a rectangular field.
#' @export

GIA <- function(v, gia.rf) {
  if (is.null(gia.rf)) {
    return(v)
  }
  n <- length(v)
  r <- gia.rf[1]
  c <- gia.rf[2]
  # corner entries
  cnpos <- c(1,c,(r-1)*c+1,r*c)
  corner <- v[cnpos]
  # geometrically isomorphic permutation
  if (which(corner == min(corner)) == 4) {
    v <- v[sort(seq(n), decreasing = TRUE)]
  }
  if (which(corner == min(corner)) == 3) {
    v <- matrix(v, r, c, byrow = TRUE)
    v <- as.vector(t(v[sort(seq(r), decreasing = TRUE), ]))
  }
  if (which(corner == min(corner)) == 2) {
    v <- matrix(v, r, c, byrow = TRUE)
    v <- as.vector(t(v[, sort(seq(c), decreasing = TRUE)]))
  }
  return(v)
}
