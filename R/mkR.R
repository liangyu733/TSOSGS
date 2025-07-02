#' @title mkR
#'
#' @description
#' A function to generate spatial correlation matrix R for a given retangular field.
#'
#' @param gia.rf A integer vector of length 2 indicating nrow and ncol of a retangular field for geometrically isomorphic allocation. See function `GIA`.
#' @param rho A numeric value for spatial correlation parameter ranges from 0 to 1.
#' @param fun 3 spatial correlation functions "power", "exponential", and "Gaussian" are provided.
#'
#' @export

mkR <- function(gia.rf, rho = 0.5, fun = "power") {
  if (!fun %in% c("power", "exponential", "Gaussian")) {
    fun == "power"
    warning("Invalid fun type; Default fun 'power' was applied.")
  }
  if (!is.vector(gia.rf) || length(gia.rf) != 2) {
    stop("gia.rf <- c(r, c)")
  }
  r <- gia.rf[1]
  c <- gia.rf[2]
  rownum <- rep(c(1:r), each = c)
  colnum <- rep(c(1:c), r)
  plotnum <- rbind(rownum, colnum)
  R <- matrix(NA, r*c, r*c)
  for (i in seq((r*c)^2)) {
    R[i] <- spcor(plotnum[, (i-1)%%(r*c)+1],plotnum[, (i-1)%/%(r*c)+1], rho, fun)
  }
  return(R)
}

